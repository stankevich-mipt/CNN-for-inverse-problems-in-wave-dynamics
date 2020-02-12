#pragma once
#include <windows.h>
class PerformanceTimer
{
public:
	PerformanceTimer()
	{
		Init();
	}
	void Reset()
	{
		timeElasped = 0;
	}
	void Start()
	{
		startTime = GetTimeCounter();
	}
	double Stop_GetDelta()
	{

		__int64 delta = GetTimeCounter() - startTime;
		timeElasped += delta;
		return delta * frequency;
	}
	double GetAccumulatedTime() const
	{
		return timeElasped * frequency;
	}
private:
	void Init()
	{
		LARGE_INTEGER freqstruct;
		QueryPerformanceFrequency( &freqstruct );
		frequency = 1.0 / double(freqstruct.QuadPart);
		Reset();
		startTime = 0;
	}
	__int64 GetTimeCounter()
	{
		LARGE_INTEGER timestruct;
		QueryPerformanceCounter( &timestruct );
		return timestruct.QuadPart;
	}
	__int64 timeElasped;
	__int64 startTime;
	double frequency;
};
