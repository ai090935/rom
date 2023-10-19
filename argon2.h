#ifndef argon2_h
#define argon2_h
#include<cstdint>

//RFC 9106

typedef unsigned char byte;

struct array
{
	byte* data;
	uint64_t length;
};

struct const_array
{
	const byte* data;
	uint64_t length;
};

struct argon2_input
{
	const_array password;
	const_array salt;
	const_array secret;
	const_array associated_data;
};

struct argon2_option
{
	uint32_t time_cost;
	uint32_t memory_cost;
	uint32_t parallelism;
};

void argon2i(argon2_input input, argon2_option option, const_array rom, array output);
void argon2d(argon2_input input, argon2_option option, const_array rom, array output);
void argon2id(argon2_input input, argon2_option option, const_array rom, array output);

#endif