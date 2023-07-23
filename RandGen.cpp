#include "RandGen.h"

bool randNum() {
	std::random_device os_seed;
	const u32 seed = os_seed();

	engine generator(seed);
	std::uniform_int_distribution<u32> distribute(0, 1);
	return distribute(generator);
}