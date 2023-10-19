#include"argon2.h"
#include<algorithm>
#include<mutex>
#include<condition_variable>
#include<vector>
#include<future>
#include<tuple>
#include<memory>
#include<optional>

namespace
{
	template<typename word>
	word byte_to_word(const byte* p) noexcept
	{
		word n = 0;
		for (uint64_t i = 0; i < sizeof(word); i++)
			n |= static_cast<word>(p[i]) << i * 8;

		return n;
	}

	template<typename word>
	void word_to_byte(word n, byte* p) noexcept
	{
		for (uint64_t i = 0; i < sizeof(word); i++)
			p[i] = static_cast<byte>(n >> i * 8);
	}

	template<typename T>
	constexpr T rotr(T value, int shift) noexcept
	{
		return (value >> shift) | (value << (sizeof(T) * 8 - shift));
	}

	template<int bits, typename T>
	constexpr T high_bits(T value) noexcept
	{
		return value >> (sizeof(T) * 8 - bits);
	}

	template<int bits, typename T>
	constexpr T low_bits(T value) noexcept
	{
		return value & high_bits<bits>(~T(0));
	}

	//-------------------------------------------------------------------------------------------------

	template<typename T>
	class matrix
	{
	public:
		matrix(size_t row, size_t col) : col(col), data(new T[row * col]) {}

		T* operator[](size_t row) noexcept
		{
			return this->data.get() + row * this->col;
		}

	private:
		size_t col;
		std::unique_ptr<T[]> data;
	};

	class barrier
	{
	public:
		barrier(size_t expected) : arrive_count(0), expected(expected), phase(0) {}
		barrier(const barrier&) = delete;
		barrier& operator=(const barrier&) = delete;

		void arrive_and_wait() noexcept
		{
			std::unique_lock lock(this->mutex);
			this->arrive_count++;

			if (this->arrive_count < this->expected)
				this->condition_variable.wait(lock, [current_phase = this->phase, this]() { return current_phase != this->phase; });
			else
			{
				this->arrive_count = 0;
				this->phase++;

				lock.unlock();
				this->condition_variable.notify_all();
			}
		}

	private:
		size_t arrive_count;
		size_t expected;
		size_t phase;
		std::mutex mutex;
		std::condition_variable condition_variable;
	};

	//-------------------------------------------------------------------------------------------------

	enum class argon2_type
	{
		argon2d = 0,
		argon2i = 1,
		argon2id = 2,
	};

	class argon2_block
	{
	public:
		argon2_block() = default;

		argon2_block& operator^=(const argon2_block& other) noexcept
		{
			for (int i = 0; i < 128; i++)
				this->data[i] ^= other.data[i];
			return *this;
		}

		uint64_t data[128];
	};

	argon2_block operator^(const argon2_block& a, const argon2_block& b) noexcept
	{
		return argon2_block(a) ^= b;
	}

	static_assert(sizeof(argon2_block) == 1024);

	class argon2_instance
	{
	public:
		argon2_instance(argon2_option option, uint32_t rom_size, argon2_type type) noexcept : passes(option.time_cost), memory_size(option.memory_cost), lanes(option.parallelism), rom_size(rom_size), memory_blocks(4 * lanes * (memory_size / (4 * lanes))), lane_length(memory_blocks / lanes), segment_length(lane_length / 4), type(type) {}

		const uint32_t passes;
		const uint32_t memory_size;
		const uint32_t lanes;
		const uint32_t rom_size;

		const uint32_t memory_blocks;
		const uint32_t lane_length;
		const uint32_t segment_length;

		const argon2_type type;
	};

	class argon2_position
	{
	public:
		argon2_position(uint32_t pass, uint32_t lane, int slice) noexcept : pass(pass), lane(lane), slice(slice) {}

		const uint32_t pass;
		const uint32_t lane;
		const int slice;
	};

	//-------------------------------------------------------------------------------------------------

	class uint128_t
	{
	public:
		constexpr uint128_t(uint64_t n) noexcept : high(0), low(n) {}

		constexpr uint128_t& operator+=(const uint128_t& other) noexcept
		{
			this->high += other.high + (static_cast<uint64_t>(this->low + other.low) < this->low);
			this->low += other.low;

			return *this;
		}

		uint64_t high;
		uint64_t low;
	};

	//RFC 7693
	class blake2b
	{
	public:
		blake2b(int output_size = 64) noexcept : h{ 0x6a09e667f3bcc908 ^ 0x01010000 ^ static_cast<uint64_t>(output_size), 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1, 0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179 }, b_length(0), output_size(output_size), total_length(0) {}

		blake2b(const byte* input, size_t length, byte* output, int output_size = 64) noexcept : blake2b(output_size)
		{
			this->update(input, length);
			this->final(output);
		}

		void update(const byte* input, size_t length) noexcept
		{
			while (length)
			{
				if (b_length == 128)
				{
					total_length += b_length;
					compress(false);
					b_length = 0;
				}

				int outlen = (b_length + length < 128) ? length : 128 - b_length;
				std::copy(input, input + outlen, b + b_length);

				b_length += outlen;
				input += outlen;
				length -= outlen;
			}
		}

		void final(byte* output) noexcept
		{
			total_length += b_length;
			std::fill(b + b_length, b + 128, 0);
			compress(true);

			byte temp[64];
			for (int i = 0; i < 8; i++)
				word_to_byte(h[i], temp + i * sizeof(uint64_t));
			std::copy(temp, temp + output_size, output);
		}

		//blake2 rfc 3.1
		void mixing(uint64_t* v, int a, int b, int c, int d, uint64_t x, uint64_t y) noexcept
		{
			v[a] = v[a] + v[b] + x;
			v[d] = rotr(static_cast<uint64_t>(v[d] ^ v[a]), 32);
			v[c] = v[c] + v[d];
			v[b] = rotr(static_cast<uint64_t>(v[b] ^ v[c]), 24);
			v[a] = v[a] + v[b] + y;
			v[d] = rotr(static_cast<uint64_t>(v[d] ^ v[a]), 16);
			v[c] = v[c] + v[d];
			v[b] = rotr(static_cast<uint64_t>(v[b] ^ v[c]), 63);
		}

		//blake2 rfc 3.2
		void compress(bool last) noexcept
		{
			constexpr byte sigma[10][16] = {
				{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
				{ 14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3 },
				{ 11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4 },
				{ 7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8 },
				{ 9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13 },
				{ 2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9 },
				{ 12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11 },
				{ 13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10 },
				{ 6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5 },
				{ 10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0 },
			};
			constexpr uint64_t IV[8] = { 0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1, 0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179 };

			uint64_t v[16], m[16];

			std::copy(h, h + 8, v);
			std::copy(IV, IV + 8, v + 8);
			v[12] ^= total_length.low;
			v[13] ^= total_length.high;
			if (last)
				v[14] = ~v[14];

			for (int i = 0; i < 16; i++)
				m[i] = byte_to_word<uint64_t>(b + i * sizeof(uint64_t));

			for (int i = 0; i < 12; i++)
			{
				mixing(v, 0, 4, 8, 12, m[sigma[i % 10][0]], m[sigma[i % 10][1]]);
				mixing(v, 1, 5, 9, 13, m[sigma[i % 10][2]], m[sigma[i % 10][3]]);
				mixing(v, 2, 6, 10, 14, m[sigma[i % 10][4]], m[sigma[i % 10][5]]);
				mixing(v, 3, 7, 11, 15, m[sigma[i % 10][6]], m[sigma[i % 10][7]]);
				mixing(v, 0, 5, 10, 15, m[sigma[i % 10][8]], m[sigma[i % 10][9]]);
				mixing(v, 1, 6, 11, 12, m[sigma[i % 10][10]], m[sigma[i % 10][11]]);
				mixing(v, 2, 7, 8, 13, m[sigma[i % 10][12]], m[sigma[i % 10][13]]);
				mixing(v, 3, 4, 9, 14, m[sigma[i % 10][14]], m[sigma[i % 10][15]]);
			}

			for (int i = 0; i < 8; i++)
				h[i] ^= v[i] ^ v[i + 8];
		}

	private:
		uint64_t h[8];
		byte b[128];
		int b_length;
		int output_size;
		uint128_t total_length;
	};

	//rfc 3.3
	void H_prime(const byte* in, uint32_t len, byte* out, uint32_t outlen) noexcept
	{
		byte temp[4];
		blake2b h(outlen <= 64 ? outlen : 64);

		word_to_byte(outlen, temp);
		h.update(temp, 4);
		h.update(in, len);
		h.final(out);

		if (outlen > 64)
		{
			byte v[64];
			const uint32_t r = outlen / 32 + static_cast<bool>(outlen % 32) - 2;

			std::copy(out, out + 64, v);
			for (uint32_t i = 1; i < r; i++)
			{
				blake2b(v, 64, v);
				std::copy(v, v + 32, out + i * 32);
			}
			blake2b(v, 64, out + r * 32, outlen - 32 * r);
		}
	}

	//-------------------------------------------------------------------------------------------------

	//rfc 3.6
	void GB(uint64_t& a, uint64_t& b, uint64_t& c, uint64_t& d) noexcept
	{
		a = a + b + 2 * low_bits<32>(a) * low_bits<32>(b);
		d = rotr(static_cast<uint64_t>(d ^ a), 32);
		c = c + d + 2 * low_bits<32>(c) * low_bits<32>(d);
		b = rotr(static_cast<uint64_t>(b ^ c), 24);

		a = a + b + 2 * low_bits<32>(a) * low_bits<32>(b);
		d = rotr(static_cast<uint64_t>(d ^ a), 16);
		c = c + d + 2 * low_bits<32>(c) * low_bits<32>(d);
		b = rotr(static_cast<uint64_t>(b ^ c), 63);
	}

	//rfc 3.6
	void P(uint64_t* s0, uint64_t* s1, uint64_t* s2, uint64_t* s3, uint64_t* s4, uint64_t* s5, uint64_t* s6, uint64_t* s7) noexcept
	{
		GB(s0[0], s2[0], s4[0], s6[0]);
		GB(s0[1], s2[1], s4[1], s6[1]);
		GB(s1[0], s3[0], s5[0], s7[0]);
		GB(s1[1], s3[1], s5[1], s7[1]);

		GB(s0[0], s2[1], s5[0], s7[1]);
		GB(s0[1], s3[0], s5[1], s6[0]);
		GB(s1[0], s3[1], s4[0], s6[1]);
		GB(s1[1], s2[0], s4[1], s7[0]);
	}

	//rfc 3.5
	argon2_block G(const argon2_block& x, const argon2_block& y) noexcept
	{
		const argon2_block r = x ^ y;

		argon2_block z = r;
		for (int i = 0; i < 8; i++)
		{
			uint64_t* row = z.data + i * 16;
			P(row + 0, row + 2, row + 4, row + 6, row + 8, row + 10, row + 12, row + 14);
		}
		for (int i = 0; i < 8; i++)
		{
			uint64_t* col = z.data + i * 2;
			P(0 * 16 + col, 1 * 16 + col, 2 * 16 + col, 3 * 16 + col, 4 * 16 + col, 5 * 16 + col, 6 * 16 + col, 7 * 16 + col);
		}

		return z ^ r;
	}

	//-------------------------------------------------------------------------------------------------

	//rfc 3.4.2
	std::tuple<uint32_t, uint32_t, std::optional<uint32_t>> mapping_index(uint32_t j1, uint32_t j2, argon2_instance instance, argon2_position position, uint32_t index) noexcept
	{
		auto mapping = [&](uint32_t w)
		{
			const uint64_t x = (static_cast<uint64_t>(j1) * static_cast<uint64_t>(j1)) >> 32;
			const uint64_t y = (w * x) >> 32;
			const uint64_t zz = w - 1 - y;

			return zz;
		};

		if (instance.rom_size && (index % 16 == 0))
			return { 0, 0, mapping(instance.rom_size) };

		const uint32_t l = (position.pass == 0 && position.slice == 0) ? position.lane : j2 % instance.lanes;

		const uint32_t finished_blocks = (position.pass == 0 ? position.slice : 3) * instance.segment_length;
		const uint32_t w = l == position.lane ? finished_blocks + index - 1 : finished_blocks - (index == 0 ? 1 : 0);

		const uint32_t start_position = (position.pass == 0 || position.slice == 3) ? 0 : (position.slice + 1) * instance.segment_length;
		const uint32_t z = (start_position + mapping(w)) % instance.lane_length;

		return { l, z, std::nullopt };
	}

	//rfc 3.4.1.1
	auto compute_argon2d_index(const argon2_block& block, argon2_instance instance, argon2_position position, uint32_t index) noexcept
	{
		uint32_t j1 = low_bits<32>(block.data[0]);
		uint32_t j2 = high_bits<32>(block.data[0]);
		return mapping_index(j1, j2, instance, position, index);
	}

	//rfc 3.4.1.2
	auto compute_argon2i_index(const argon2_block& block, argon2_instance instance, argon2_position position, uint32_t index) noexcept
	{
		uint32_t j1 = low_bits<32>(block.data[index % 128]);
		uint32_t j2 = high_bits<32>(block.data[index % 128]);
		return mapping_index(j1, j2, instance, position, index);
	}

	void init_argon2i_input_block(argon2_block& input_block, argon2_instance instance, argon2_position position) noexcept
	{
		input_block.data[0] = position.pass;
		input_block.data[1] = position.lane;
		input_block.data[2] = position.slice;
		input_block.data[3] = instance.memory_blocks;
		input_block.data[4] = instance.passes;
		input_block.data[5] = static_cast<uint64_t>(instance.type);
		std::fill(input_block.data + 6, input_block.data + 128, 0);
	}

	void compute_argon2i_index_block(argon2_block& input_block, argon2_block& index_block) noexcept
	{
		constexpr argon2_block zero_block = {};

		input_block.data[6]++;
		index_block = G(zero_block, G(zero_block, input_block));
	}

	void compute_segment(matrix<argon2_block>& B, const argon2_block* rom, argon2_instance instance, argon2_position position) noexcept
	{
		argon2_block input_block, index_block;
		const bool is_argon2i_index = (instance.type == argon2_type::argon2i) || (instance.type == argon2_type::argon2id && position.pass == 0 && position.slice < 2);

		if (is_argon2i_index)
		{
			init_argon2i_input_block(input_block, instance, position);
			compute_argon2i_index_block(input_block, index_block);
		}

		for (uint32_t index = (position.pass == 0 && position.slice == 0 ? 2 : 0); index < instance.segment_length; index++)
		{
			const uint32_t absolute_index = instance.segment_length * position.slice + index;
			argon2_block& curr = B[position.lane][absolute_index];
			argon2_block& prev = B[position.lane][absolute_index != 0 ? absolute_index - 1 : instance.lane_length - 1];

			if (is_argon2i_index && index % 128 == 0 && index != 0)
				compute_argon2i_index_block(input_block, index_block);

			auto [l, z, rom_index] = is_argon2i_index ? compute_argon2i_index(index_block, instance, position, index) : compute_argon2d_index(prev, instance, position, index);
			const argon2_block& ref = !rom_index.has_value() ? B[l][z] : rom[rom_index.value()];

			if (position.pass == 0)
				curr = G(prev, ref);
			else
				curr ^= G(prev, ref);
		}
	}

	//-------------------------------------------------------------------------------------------------

	//rfc 3.2.1
	void compute_H0(argon2_input input, argon2_option option, array output, argon2_type type, byte* H0) noexcept
	{
		blake2b h;
		byte temp[sizeof(uint32_t)];

		word_to_byte(option.parallelism, temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(static_cast<uint32_t>(output.length), temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(option.memory_cost, temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(option.time_cost, temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(static_cast<uint32_t>(0x13), temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(static_cast<uint32_t>(type), temp);
		h.update(temp, sizeof(uint32_t));

		word_to_byte(static_cast<uint32_t>(input.password.length), temp);
		h.update(temp, sizeof(uint32_t));

		h.update(input.password.data, input.password.length);

		word_to_byte(static_cast<uint32_t>(input.salt.length), temp);
		h.update(temp, sizeof(uint32_t));

		h.update(input.salt.data, input.salt.length);

		word_to_byte(static_cast<uint32_t>(input.secret.length), temp);
		h.update(temp, sizeof(uint32_t));

		h.update(input.secret.data, input.secret.length);

		word_to_byte(static_cast<uint32_t>(input.associated_data.length), temp);
		h.update(temp, sizeof(uint32_t));

		h.update(input.associated_data.data, input.associated_data.length);

		h.final(H0);
	}

	//rfc 3.2.1 ~ 3.2.4
	void init(matrix<argon2_block>& B, argon2_input input, argon2_option option, array output, argon2_type type) noexcept
	{
		byte in[72], out[1024];

		compute_H0(input, option, output, type, in);
		std::fill(in + 64, in + 68, 0);
		for (uint32_t i = 0; i < option.parallelism; i++)
		{
			word_to_byte(i, in + 68);
			H_prime(in, 72, out, 1024);
			for (size_t j = 0; j < 128; j++)
				B[i][0].data[j] = byte_to_word<uint64_t>(out + j * sizeof(uint64_t));
		}

		in[64]++;
		for (uint32_t i = 0; i < option.parallelism; i++)
		{
			word_to_byte(i, in + 68);
			H_prime(in, 72, out, 1024);
			for (size_t j = 0; j < 128; j++)
				B[i][1].data[j] = byte_to_word<uint64_t>(out + j * sizeof(uint64_t));
		}
	}

	//rfc 3.2.5 ~ 3.2.6
	void iterator(matrix<argon2_block>& B, const argon2_block* rom, argon2_instance instance)
	{
		bool stop = true;
		std::mutex mutex;
		barrier sync_point(instance.lanes);
		std::vector<std::future<void>> vec(instance.lanes);

		auto f = [&](uint32_t lane)
		{
			mutex.lock();
			mutex.unlock();
			if (stop)
				return;

			for (uint32_t pass = 0; pass < instance.passes; pass++)
				for (int slice = 0; slice < 4; slice++)
				{
					compute_segment(B, rom, instance, argon2_position(pass, lane, slice));
					sync_point.arrive_and_wait();
				}
		};

		try
		{
			mutex.lock();
			for (uint32_t lane = 0; lane < instance.lanes; lane++)
				vec[lane] = std::async(std::launch::async, f, lane);
		}
		catch (...)
		{
			mutex.unlock();
			for (const auto& i : vec)
				if (i.valid())
					i.wait();
			throw;
		}

		stop = false;
		mutex.unlock();
		for (const auto& i : vec)
			i.wait();
	}

	//rfc 3.2.7
	void final(matrix<argon2_block>& B, uint32_t row, uint32_t col, array output) noexcept
	{
		argon2_block& C = B[0][col - 1];
		for (uint32_t i = 1; i < row; i++)
			C ^= B[i][col - 1];

		byte temp[1024];
		for (size_t i = 0; i < 128; i++)
			word_to_byte(C.data[i], temp + i * sizeof(uint64_t));
		H_prime(temp, 1024, output.data, output.length);
	}

	void argon2(argon2_input input, argon2_option option, const_array rom, array output, argon2_type type)
	{
		argon2_instance instance(option, rom.length / 1024, type);
		matrix<argon2_block> B(instance.lanes, instance.lane_length);

		init(B, input, option, output, type);
		iterator(B, reinterpret_cast<const argon2_block*>(rom.data), instance);
		final(B, instance.lanes, instance.lane_length, output);
	}
}

//-------------------------------------------------------------------------------------------------

void argon2i(argon2_input input, argon2_option option, const_array rom, array output)
{
	argon2(input, option, rom, output, argon2_type::argon2i);
}

void argon2d(argon2_input input, argon2_option option, const_array rom, array output)
{
	argon2(input, option, rom, output, argon2_type::argon2d);
}

void argon2id(argon2_input input, argon2_option option, const_array rom, array output)
{
	argon2(input, option, rom, output, argon2_type::argon2id);
}