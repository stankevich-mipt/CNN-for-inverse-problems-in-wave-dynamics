#ifndef AleaneraRelitningTest3D120926
#define AleaneraRelitningTest3D120926

#include <stdlib.h>


void Test3D()
//2 тетраэдра. —торона 012 обща€, вершины 3,4 принадлежат различным.
{
  RealType * R[5];
  RealType V_n[5];//node
  RealType * V_e[9];//edge
  RealType * V_s[7];//side
  RealType * V_v[2];//volume
  RealType x, y, z;
  int a, b, c, A, B, C, D;
  int i, j, k;
  int N;
  int N_e, N_s, N_v;
  RealType * Ri;
  RealType * Vi;
  RealType Result[3];
  char I0, I1;
  int Num;
  RealType Error;
  RealType Inter;
  RealType Func;

  scanf("%i",&Num);

  //инициализируем координаты точек
  for(i=0; (i<5); i++)
  {
    R[i]=malloc(3*sizeof(RealType));
    if((R[i])==NULL)
    {
      printf("Malloc or realloc error\n");
      ErrorExit();
    };
  };   
  (R[0])[0]=0;
  (R[0])[1]=0;
  (R[0])[2]=0;

  (R[1])[0]=1;
  (R[1])[1]=1;
  (R[1])[2]=0;

  (R[2])[0]=0;
  (R[2])[1]=0;
  (R[2])[2]=1;

  (R[3])[0]=0;
  (R[3])[1]=1;
  (R[3])[2]=0;

  (R[4])[0]=1;
  (R[4])[1]=0;
  (R[4])[2]=0;

  I0='F';
  if((interpolationCompositionalProduct3D(
(R[2])[0]-(R[3])[0], (R[2])[1]-(R[3])[1], (R[2])[2]-(R[3])[2],
(R[0])[0]-(R[3])[0], (R[0])[1]-(R[3])[1], (R[0])[2]-(R[3])[2],
(R[1])[0]-(R[3])[0], (R[1])[1]-(R[3])[1], (R[1])[2]-(R[3])[2]))>0)
  {
    I0='T';
  };

  I1='F';
  if((interpolationCompositionalProduct3D(
(R[2])[0]-(R[4])[0], (R[2])[1]-(R[4])[1], (R[2])[2]-(R[4])[2],
(R[0])[0]-(R[4])[0], (R[0])[1]-(R[4])[1], (R[0])[2]-(R[4])[2],
(R[1])[0]-(R[4])[0], (R[1])[1]-(R[4])[1], (R[1])[2]-(R[4])[2]))>0)
  {
    I1='T';
  };

//цикл тестировани€ по N
  for(N=1; (N<6); N++)
  {
    N_e=N-1;
    N_s=((N-1)*(N-2))/2;
    N_v=((N-1)*(N-2)*(N-3))/6;
    if(N_e>0)
    {
      for(i=0; (i<9); i++)
      {
        V_e[i]=malloc(N_e*sizeof(RealType));
        if((V_e[i])==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
      };
    };
    if(N_s>0)
    {
      for(i=0; (i<7); i++)
      {
        V_s[i]=malloc(N_s*sizeof(RealType));
        if((V_s[i])==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
      };
    };
    if(N_v>0)
    {
      for(i=0; (i<2); i++)
      {
        V_v[i]=malloc(N_v*sizeof(RealType));
        if((V_v[i])==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
      };
    };
//цикл по степени полинома
    for(GMFunc=0; (GMFunc<=(N+1)); GMFunc++)
    {
//инициализируем опорные точки значени€ми полинома
      for(i=0; (i<5); i++)
      {
        V_n[i]=func_3D((R[i])[0], (R[i])[1], (R[i])[2]);
      };

      if(N_e>0)
      {
        Ri=malloc(6*sizeof(RealType));
        if(Ri==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
        for(i=0; (i<9); i++)
        {
          if(i==0) {A=0; B=1;}
          else if(i==1) {A=0; B=2;}
          else if(i==2) {A=0; B=3;}
          else if(i==3) {A=0; B=4;}
          else if(i==4) {A=1; B=2;}
          else if(i==5) {A=1; B=3;}
          else if(i==6) {A=1; B=4;}
          else if(i==7) {A=2; B=3;}
          else {A=2; B=4;};
          for(j=0; (j<3); j++)
          {
            Ri[j]=(R[A])[j];
            Ri[3+j]=(R[B])[j];
          };
          for(j=0; (j<N_e); j++)
          {
            Interpolation_3D_Get_Edge(Result, N, Ri, j);
            (V_e[i])[j]=func_3D(Result[0], Result[1], Result[2]);
          };
        };
        free(Ri);
      };

      if(N_s>0)
      {
        Ri=malloc(9*sizeof(RealType));
        if(Ri==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
        for(i=0; (i<7); i++)
        {
          if(i==0) {A=0; B=1; C=2;}
          else if(i==1) {A=0; B=1; C=3;}
          else if(i==2) {A=0; B=1; C=4;}
          else if(i==3) {A=0; B=2; C=3;}
          else if(i==4) {A=0; B=2; C=4;}
          else if(i==5) {A=1; B=2; C=3;}
          else {A=1; B=2; C=4;};
          for(j=0; (j<3); j++)
          {
            Ri[j]=(R[A])[j];
            Ri[3+j]=(R[B])[j];
            Ri[6+j]=(R[C])[j];
          };
          for(j=0; (j<N_s); j++)
          {
            Interpolation_3D_Get_Side(Result, N, Ri, j);
            (V_s[i])[j]=func_3D(Result[0], Result[1], Result[2]);
          };
        };
        free(Ri);
      };
        
      if(N_v>0)
      {
        Ri=malloc(12*sizeof(RealType));
        if(Ri==NULL)
        {
          printf("Malloc or realloc error\n");
          ErrorExit();
        };
        for(i=0; (i<2); i++)
        {
          if(i==0) {A=0; B=1; C=2; D=3;}
          else {A=0; B=1; C=2; D=4;};
          for(j=0; (j<3); j++)
          {
            Ri[j]=(R[A])[j];
            Ri[3+j]=(R[B])[j];
            Ri[6+j]=(R[C])[j];
            Ri[9+j]=(R[D])[j];
          };
          for(j=0; (j<N_v); j++)
          {
            Interpolation_3D_Get_Volume(Result, N, Ri, j);
            (V_v[i])[j]=func_3D(Result[0], Result[1], Result[2]);
          };
        };
        free(Ri);
      };

//цикл по специальным тестовым опорным точкам
      Error=0;
      Ri=malloc(15*sizeof(RealType));
      if(Ri==NULL)
      {
        printf("Malloc or realloc error\n");
        ErrorExit();
      };
      Vi=malloc((((N+1)*(N+2)*(N+3))/6)*sizeof(RealType));
      if(Vi==NULL)
      {
        printf("Malloc or realloc error\n");
        ErrorExit();
      };
      for(a=0; (a<=Num); a++)
      {
        for(b=0; ((a+b)<=Num); b++)
        {
          for(c=0; ((a+b+c)<=Num); c++)
          {
//первый тетраэдр
            A=0; B=1; C=2; D=3;
            x=(R[D])[0]+((((RealType)a)/((RealType)Num))*((R[A])[0]-(R[D])[0])+(((RealType)b)/((RealType)Num))*((R[B])[0]-(R[D])[0])+(((RealType)c)/((RealType)Num))*((R[C])[0]-(R[D])[0]));
            y=(R[D])[1]+((((RealType)a)/((RealType)Num))*((R[A])[1]-(R[D])[1])+(((RealType)b)/((RealType)Num))*((R[B])[1]-(R[D])[1])+(((RealType)c)/((RealType)Num))*((R[C])[1]-(R[D])[1]));
            z=(R[D])[2]+((((RealType)a)/((RealType)Num))*((R[A])[2]-(R[D])[2])+(((RealType)b)/((RealType)Num))*((R[B])[2]-(R[D])[2])+(((RealType)c)/((RealType)Num))*((R[C])[2]-(R[D])[2]));
            Func=func_3D(x, y, z);
            Ri[0]=x;
            Ri[1]=y;
            Ri[2]=z;
            for(i=0; (i<3); i++)
            {
              Ri[3+i]=(R[A])[i];
              Ri[6+i]=(R[B])[i];
              Ri[9+i]=(R[C])[i];
              Ri[12+i]=(R[D])[i];
            };
            Vi[0]=V_n[A];
            Vi[1]=V_n[B];
            Vi[2]=V_n[C];
            Vi[3]=V_n[D];
            k=4;
            if(N_e>0)
            {
              for(i=0; (i<9); i++)
              {
                if((i==0)||(i==1)||(i==2)||(i==4)||(i==5)||(i==7))
                {
                  for(j=0; (j<N_e); j++, k++)
                  {
                    Vi[k]=(V_e[i])[j];
                  };
                };
              };
            };
            if(N_s>0)
            {
              for(i=0; (i<7); i++)
              {
                if((i==0)||(i==1)||(i==3)||(i==5))
                {
                  for(j=0; (j<N_s); j++, k++)
                  {
                    Vi[k]=(V_s[i])[j];
                  };
                };
              };
            };
            if(N_v>0)
            {
              i=0;
              for(j=0; (j<N_v); j++, k++)
              {
                Vi[k]=(V_v[i])[j];
              };
            };
            Inter=Interpolate_3D(N, 'P', Ri, I0, Vi);
            Error=(Inter-Func)*(Inter-Func)+Error;
//второй тетраэдр
            A=0; B=1; C=2; D=4;
            x=(R[D])[0]+((((RealType)a)/((RealType)Num))*((R[A])[0]-(R[D])[0])+(((RealType)b)/((RealType)Num))*((R[B])[0]-(R[D])[0])+(((RealType)c)/((RealType)Num))*((R[C])[0]-(R[D])[0]));
            y=(R[D])[1]+((((RealType)a)/((RealType)Num))*((R[A])[1]-(R[D])[1])+(((RealType)b)/((RealType)Num))*((R[B])[1]-(R[D])[1])+(((RealType)c)/((RealType)Num))*((R[C])[1]-(R[D])[1]));
            z=(R[D])[2]+((((RealType)a)/((RealType)Num))*((R[A])[2]-(R[D])[2])+(((RealType)b)/((RealType)Num))*((R[B])[2]-(R[D])[2])+(((RealType)c)/((RealType)Num))*((R[C])[2]-(R[D])[2]));
            Func=func_3D(x, y, z);
            Ri[0]=x;
            Ri[1]=y;
            Ri[2]=z;
            for(i=0; (i<3); i++)
            {
              Ri[3+i]=(R[A])[i];
              Ri[6+i]=(R[B])[i];
              Ri[9+i]=(R[C])[i];
              Ri[12+i]=(R[D])[i];
            };
            Vi[0]=V_n[A];
            Vi[1]=V_n[B];
            Vi[2]=V_n[C];
            Vi[3]=V_n[D];
            k=4;
            if(N_e>0)
            {
              for(i=0; (i<9); i++)
              {
                if((i==0)||(i==1)||(i==3)||(i==4)||(i==6)||(i==8))
                {
                  for(j=0; (j<N_e); j++, k++)
                  {
                    Vi[k]=(V_e[i])[j];
                  };
                };
              };
            };
            if(N_s>0)
            {
              for(i=0; (i<7); i++)
              {
                if((i==0)||(i==2)||(i==4)||(i==6))
                {
                  for(j=0; (j<N_s); j++, k++)
                  {
                    Vi[k]=(V_s[i])[j];
                  };
                };
              };
            };
            if(N_v>0)
            {
              i=1;
              for(j=0; (j<N_v); j++, k++)
              {
                Vi[k]=(V_v[i])[j];
              };
            };
            Inter=Interpolate_3D(N, 'P', Ri, I1, Vi);
            Error=(Inter-Func)*(Inter-Func)+Error;
          };
        };
      };
      free(Ri);
      free(Vi);
      Error=Error/((RealType)(((Num+1)*(Num+2)*(Num+3))/6));
      printf("N=%i F=%i Error=%f\n", N, GMFunc, Error);
    };
    if(N_e>0)
    {
      for(i=0; (i<9); i++)
      {
        free(V_e[i]);
      };
    };
    if(N_s>0)
    {
      for(i=0; (i<7); i++)
      {
        free(V_s[i]);
      };
    };
    if(N_v>0)
    {
      for(i=0; (i<2); i++)
      {
        free(V_v[i]);
      };
    };
  };
  
  for(i=0; (i<5); i++)
  {
    free(R[i]);
  };   
  return;
}

#endif


