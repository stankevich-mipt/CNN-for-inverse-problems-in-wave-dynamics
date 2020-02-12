// #define _HAS_EXCEPTIONS 0
// #define _ITERATOR_DEBUG_LEVEL 0 
// #include <xmmintrin.h>
// #include <vld.h>

#include "Task.h"

//#define ORDER_FROM_SETTINGS
#define SPACE_FROM_SETTINGS

const unsigned int defaultPolynomialOrder = 1;
typedef Space2 DefaultSpace;

int main()
{
  #if(defined(SPACE_FROM_SETTINGS) && defined(ORDER_FROM_SETTINGS))
  {
    BasicSettings basicSettings; 
    basicSettings.Parse("task.xml");

    switch(basicSettings.configDimsCount)
    {
      case 2:
      {
        Settings<Space2> settings;
        settings.Parse("task.xml");

        switch (settings.solver.polynomialsOrder) 
        {
          case 1: {Task<Space2, 1> task; task.Run();} break;
          case 2: {Task<Space2, 2> task; task.Run();} break;
          case 3: {Task<Space2, 3> task; task.Run();} break;
          case 4: {Task<Space2, 4> task; task.Run();} break;
          case 5: {Task<Space2, 5> task; task.Run();} break;
          case 6: {Task<Space2, 6> task; task.Run();} break;
          case 7: {Task<Space2, 7> task; task.Run();} break;
          //case 8: {Task<Space2, 8> task; task.Run();} break;
          default: std::cerr << "Unknown polynomial order"; break;
        }
      }break;
      case 3:
      {
        Settings<Space3> settings;
        settings.Parse("task.xml");

        switch (settings.solver.polynomialsOrder) 
        {
          case 1: {Task<Space3, 1> task; task.Run();} break;
          case 2: {Task<Space3, 2> task; task.Run();} break;
          case 3: {Task<Space3, 3> task; task.Run();} break;
          case 4: {Task<Space3, 4> task; task.Run();} break;
          case 5: {Task<Space3, 5> task; task.Run();} break;
          case 6: {Task<Space3, 6> task; task.Run();} break;
          case 7: {Task<Space3, 7> task; task.Run();} break;
          //case 8: {Task<Space3, 8> task; task.Run();} break;
          default: std::cerr << "Unknown polynomial order"; break;
        }
      }break;
      default: std::cerr << "Unknown dims count"; break;
    }
  }
  #endif
  #if(defined(SPACE_FROM_SETTINGS) && !defined(ORDER_FROM_SETTINGS))
  {
    BasicSettings settings; 
    settings.Parse("task.xml");

    switch(settings.configDimsCount)
    {
      case 2:
      {
        Task<Space2, defaultPolynomialOrder> task;
        task.Run(); 
      }break;
      case 3:
      {
        Task<Space3, defaultPolynomialOrder> task;
        task.Run(); 
      }break;
      default: std::cerr << "Unknown dims count"; break;
    }
  }
  #endif
  #if(!defined(SPACE_FROM_SETTINGS) && defined(ORDER_FROM_SETTINGS))
  {
    BasicSettings settings; 
    settings.Parse("task.xml");

    switch(settings.configDimsCount)
    {
      case 2:
      {
        switch (settings.solver.polynomialsOrder) 
        {
          case 1: {Task<DefaultSpace, 1> task; task.Run();} break;
          case 2: {Task<DefaultSpace, 2> task; task.Run();} break;
          case 3: {Task<DefaultSpace, 3> task; task.Run();} break;
          case 4: {Task<DefaultSpace, 4> task; task.Run();} break;
          case 5: {Task<DefaultSpace, 5> task; task.Run();} break;
          case 6: {Task<DefaultSpace, 6> task; task.Run();} break;
          case 7: {Task<DefaultSpace, 7> task; task.Run();} break;
          case 8: {Task<DefaultSpace, 8> task; task.Run();} break;
          default: std::cerr << "Unknown polynomial order"; break;
        }
      }break;
    }
  }
  #endif
  #if(!defined(SPACE_FROM_SETTINGS) && !defined(ORDER_FROM_SETTINGS))
  {
    Task<DefaultSpace, defaultPolynomialOrder> task;
    task.Run();
  }
  #endif

  return 0;
}
