#include <stdlib.h>

template<typename T>
T random(T min, T max)
{
	return min - (min - max) * (T(rand()) / T(32768));
}
