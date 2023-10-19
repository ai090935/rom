#include<iostream>
#include<chrono>
#include<fcntl.h>
#include<sys/mman.h>
#include"argon2.h"


const byte* get_rom(const char* path, size_t length)
{
	int fd = open(path, O_RDONLY);
	if (fd == -1)
		throw std::runtime_error("fd == -1");
	const void* p = mmap(nullptr, length, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
	if (p == MAP_FAILED)
		throw std::runtime_error("MAP_FAILED");

	return static_cast<const byte*>(p);
}

template<typename F>
auto measure(F&& f)
{
	auto start = std::chrono::steady_clock::now();
	f();
	auto end = std::chrono::steady_clock::now();

	return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

void print_result(const byte* p)
{
	std::cout << "result: ";
	std::cout << std::hex;
	for (int i = 0; i < 64; i++)
		std::cout << static_cast<int>(p[i] >> 4) << static_cast<int>(p[i] & 0xf);
	std::cout << std::dec << '\n';
}

template<typename F>
void run(F&& f, uint32_t t, uint32_t m, uint32_t p, const byte* rom, uint32_t rom_size)
{
	byte out[64];
	auto fn = [&]()
	{
		f({}, { t, m, p }, { rom, static_cast<uint64_t>(rom_size) * 1024 }, { out, 64 });
	};

	std::cout << "t=" << t << ", m=" << m << ", p=" << p << ", rom=" << rom_size << ": " << measure(fn) << " microseconds\n";
	print_result(out);
}

constexpr uint64_t K = 1 << 10, M = 1 << 20, G = 1 << 30;

int main()
{
	const byte* rom = get_rom("/mnt/rom", 16 * G);

	uint32_t t = 1;
	uint32_t m = 1 * M;
	uint32_t p = 4;
	run(argon2id, t, m, p, rom, 0 * M);
	//run(argon2id, t, m, p, rom, 1 * M);
	//run(argon2id, t, m, p, rom, 2 * M);
	//run(argon2id, t, m, p, rom, 4 * M);
	//run(argon2id, t, m, p, rom, 8 * M);
	//run(argon2id, t, m, p, rom, 16 * M);

	return 0;
}
