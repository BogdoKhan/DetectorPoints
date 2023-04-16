#pragma once
#include <cctype>
#include <random>

using u32 = uint_least32_t;
using engine = std::mt19937;

bool randNum() {
	std::random_device os_seed;
	const u32 seed = os_seed();

	engine generator(seed);
	std::uniform_int_distribution<u32> distribute(0, 1);
	return distribute(generator);
}