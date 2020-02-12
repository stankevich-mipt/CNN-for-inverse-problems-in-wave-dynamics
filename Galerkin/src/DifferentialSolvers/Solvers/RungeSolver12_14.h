#include "../DifferentialSolver.h"
#include <assert.h>

/*
Copy-past of DormandPrince method & light differences
results to Runge-Kutta 14th order method
http://sce.uhcl.edu/rungekutta/
http://sce.uhcl.edu/rungekutta/rk1412.txt
*/
template<typename Scalar>
class RungeSolver12_14 : public DifferentialSolver<Scalar>
{
public:
	RungeSolver12_14()
	{
	}
	void SetSystem(DifferentialSystem<Scalar> *system)
	{
		this->system = system;
		int maxDimentionsCount = system->GetMaxDimentionsCount();
		currCoords = new Scalar[maxDimentionsCount];

		nextCoords1 = new Scalar[maxDimentionsCount];
		nextCoords2 = new Scalar[maxDimentionsCount];
		derivatives = new Scalar[maxDimentionsCount];
		probeCoords = new Scalar[maxDimentionsCount];

		k1 = new Scalar[maxDimentionsCount];
		k2 = new Scalar[maxDimentionsCount];
		k3 = new Scalar[maxDimentionsCount];
		k4 = new Scalar[maxDimentionsCount];
		k5 = new Scalar[maxDimentionsCount];
		k6 = new Scalar[maxDimentionsCount];
		k7 = new Scalar[maxDimentionsCount];
		k8 = new Scalar[maxDimentionsCount];
		k9 = new Scalar[maxDimentionsCount];
		k10 = new Scalar[maxDimentionsCount];
		k11 = new Scalar[maxDimentionsCount];
		k12 = new Scalar[maxDimentionsCount];
		k13 = new Scalar[maxDimentionsCount];
		k14 = new Scalar[maxDimentionsCount];
		k15 = new Scalar[maxDimentionsCount];
		k16 = new Scalar[maxDimentionsCount];
		k17 = new Scalar[maxDimentionsCount];
		k18 = new Scalar[maxDimentionsCount];
		k19 = new Scalar[maxDimentionsCount];
		k20 = new Scalar[maxDimentionsCount];
		k21 = new Scalar[maxDimentionsCount];
		k22 = new Scalar[maxDimentionsCount];
		k23 = new Scalar[maxDimentionsCount];
		k24 = new Scalar[maxDimentionsCount];
		k25 = new Scalar[maxDimentionsCount];
		k26 = new Scalar[maxDimentionsCount];
		k27 = new Scalar[maxDimentionsCount];
		k28 = new Scalar[maxDimentionsCount];
		k29 = new Scalar[maxDimentionsCount];
		k30 = new Scalar[maxDimentionsCount];
		k31 = new Scalar[maxDimentionsCount];
		k32 = new Scalar[maxDimentionsCount];
		k33 = new Scalar[maxDimentionsCount];
		k34 = new Scalar[maxDimentionsCount];
		k35 = new Scalar[maxDimentionsCount];

		this->currTime = 0;
	}
	virtual ~RungeSolver12_14()
	{
		delete[] currCoords;
		delete[] nextCoords1;
		delete[] nextCoords2;
		delete[] derivatives;
		delete[] probeCoords;

		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		delete[] k5;
		delete[] k6;
		delete[] k7;
		delete[] k8;
		delete[] k9;
		delete[] k10;
		delete[] k11;
		delete[] k12;
		delete[] k13;
		delete[] k14;
		delete[] k15;
		delete[] k16;
		delete[] k17;
		delete[] k18;
		delete[] k19;
		delete[] k20;
		delete[] k21;
		delete[] k22;
		delete[] k23;
		delete[] k24;
		delete[] k25;
		delete[] k26;
		delete[] k27;
		delete[] k28;
		delete[] k29;
		delete[] k30;
		delete[] k31;
		delete[] k32;
		delete[] k33;
		delete[] k34;
		delete[] k35;
	}

	Scalar pow(Scalar a, Scalar b)
	{
		return exp(log(a) * b);
	}

	int GetPhasesCount()
	{
		return 35;
	}

	void InitStep(Scalar timeStep, Scalar tolerance,
		int globalStepIndex, int hierarchyPhase)
	{
		assert(system->GetHierarchyLevelsCount() == 1);
		DifferentialSolver<Scalar>::InitStep(timeStep, tolerance, globalStepIndex, hierarchyPhase);
		system->GetCurrCoords(this->currTime, currCoords, globalStepIndex, hierarchyPhase);
	}

	virtual void AdvancePhase(int phaseIndex)
	{
		switch (phaseIndex)
		{
		case 0:
		{
				  Scalar a0(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b1_0(0.111111111111111111111111111111111111111111111111111111111111);
				  system->GetCurrDerivatives(k1, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b1_0 * k1[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a0, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 1:
		{
				  Scalar a1(0.111111111111111111111111111111111111111111111111111111111111);
				  Scalar b2_0(-0.833333333333333333333333333333333333333333333333333333333333);
				  Scalar b2_1(1.38888888888888888888888888888888888888888888888888888888889);
				  system->GetCurrDerivatives(k2, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b2_0 * k1[coordIndex]
						  + b2_1 * k2[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a1, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 2:
		{
				  Scalar a2(0.555555555555555555555555555555555555555555555555555555555556);
				  Scalar b3_0(0.208333333333333333333333333333333333333333333333333333333333);
				  Scalar b3_1(0);
				  Scalar b3_2(0.625000000000000000000000000000000000000000000000000000000000);
				  system->GetCurrDerivatives(k3, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b3_0 * k1[coordIndex]
						  + b3_1 * k2[coordIndex]
						  + b3_2 * k3[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a2, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 3:
		{
				  Scalar a3(0.833333333333333333333333333333333333333333333333333333333333);
				  Scalar b4_0(0.193333333333333333333333333333333333333333333333333333333333);
				  Scalar b4_1(0);
				  Scalar b4_2(0.220000000000000000000000000000000000000000000000000000000000);
				  Scalar b4_3(-0.0800000000000000000000000000000000000000000000000000000000000);
				  system->GetCurrDerivatives(k4, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b4_0 * k1[coordIndex]
						  + b4_1 * k2[coordIndex]
						  + b4_2 * k3[coordIndex]
						  + b4_3 * k4[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a3, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 4:
		{
				  Scalar a4(0.333333333333333333333333333333333333333333333333333333333333);
				  Scalar b5_0(0.100000000000000000000000000000000000000000000000000000000000);
				  Scalar b5_1(0);
				  Scalar b5_2(0);
				  Scalar b5_3(0.400000000000000000000000000000000000000000000000000000000000);
				  Scalar b5_4(0.500000000000000000000000000000000000000000000000000000000000);
				  system->GetCurrDerivatives(k5, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b5_0 * k1[coordIndex]
						  + b5_1 * k2[coordIndex]
						  + b5_2 * k3[coordIndex]
						  + b5_3 * k4[coordIndex]
						  + b5_4 * k5[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a4, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 5:
		{
				  Scalar a5(1.00000000000000000000000000000000000000000000000000000000000);
				  Scalar b6_0(0.103484561636679776672993546511910344499744798201971316606663);
				  Scalar b6_1(0);
				  Scalar b6_2(0);
				  Scalar b6_3(0.122068887306407222589644082868962077139592714834162134741275);
				  Scalar b6_4(0.482574490331246622475134780125688112865919023850168049679402);
				  Scalar b6_5(-0.0381409600015606999730886240005620205664113072478411477421970);
				  system->GetCurrDerivatives(k6, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b6_0 * k1[coordIndex]
						  + b6_1 * k2[coordIndex]
						  + b6_2 * k3[coordIndex]
						  + b6_3 * k4[coordIndex]
						  + b6_4 * k5[coordIndex]
						  + b6_5 * k6[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a5, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 6:
		{
				  Scalar a6(0.669986979272772921764683785505998513938845229638460353285142);
				  Scalar b7_0(0.124380526654094412881516420868799316268491466359671423163289);
				  Scalar b7_1(0);
				  Scalar b7_2(0);
				  Scalar b7_3(0);
				  Scalar b7_4(0.226120282197584301422238662979202901196752320742633143965145);
				  Scalar b7_5(0.0137885887618080880607695837016477814530969417491493385363543);
				  Scalar b7_6(-0.0672210133996684449749399507414305856950086341525382182856200);
				  system->GetCurrDerivatives(k7, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b7_0 * k1[coordIndex]
						  + b7_1 * k2[coordIndex]
						  + b7_2 * k3[coordIndex]
						  + b7_3 * k4[coordIndex]
						  + b7_4 * k5[coordIndex]
						  + b7_5 * k6[coordIndex]
						  + b7_6 * k7[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a6, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 7:
		{
				  Scalar a7(0.297068384213818357389584716808219413223332094698915687379168);
				  Scalar b8_0(0.0936919065659673815530885456083005933866349695217750085655603);
				  Scalar b8_1(0);
				  Scalar b8_2(0);
				  Scalar b8_3(0);
				  Scalar b8_4(0);
				  Scalar b8_5(-0.00613406843450510987229498995641664735620914507128858871007099);
				  Scalar b8_6(0.216019825625503063708860097659866573490979433278117320188668);
				  Scalar b8_7(0.423695063515761937337619073960976753205867469544123532683116);
				  system->GetCurrDerivatives(k8, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b8_0 * k1[coordIndex]
						  + b8_1 * k2[coordIndex]
						  + b8_2 * k3[coordIndex]
						  + b8_3 * k4[coordIndex]
						  + b8_4 * k5[coordIndex]
						  + b8_5 * k6[coordIndex]
						  + b8_6 * k7[coordIndex]
						  + b8_7 * k8[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a7, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 8:
		{
				  Scalar a8(0.727272727272727272727272727272727272727272727272727272727273);
				  Scalar b9_0(0.0838479812409052664616968791372814085980533139224911131069335);
				  Scalar b9_1(0);
				  Scalar b9_2(0);
				  Scalar b9_3(0);
				  Scalar b9_4(0);
				  Scalar b9_5(-0.0117949367100973814319755056031295775367961960590736150777613);
				  Scalar b9_6(-0.247299020568812652339473838743194598325992840353340132697498);
				  Scalar b9_7(0.0978080858367729012259313014081291665503740655476733940756599);
				  Scalar b9_8(0.217590689243420631360008651767860318344168120024782176879989);
				  system->GetCurrDerivatives(k9, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b9_0 * k1[coordIndex]
						  + b9_1 * k2[coordIndex]
						  + b9_2 * k3[coordIndex]
						  + b9_3 * k4[coordIndex]
						  + b9_4 * k5[coordIndex]
						  + b9_5 * k6[coordIndex]
						  + b9_6 * k7[coordIndex]
						  + b9_7 * k8[coordIndex]
						  + b9_8 * k9[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a8, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 9:
		{
				  Scalar a9(0.140152799042188765276187487966946717629806463082532936287323);
				  Scalar b10_0(0.0615255359769428227954562389614314714333423969064821107453940);
				  Scalar b10_1(0);
				  Scalar b10_2(0);
				  Scalar b10_3(0);
				  Scalar b10_4(0);
				  Scalar b10_5(0.00592232780324503308042990005798046524738389560444257136834990);
				  Scalar b10_6(0.470326159963841112217224303205894113455362530746108825010848);
				  Scalar b10_7(0.299688863848679000853981837096192399136831121671781279184194);
				  Scalar b10_8(-0.247656877593994914689992276329810825853958069263947095548189);
				  Scalar b10_9(0.110895029771437682893999851839061714522445173600678718208625);
				  system->GetCurrDerivatives(k10, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b10_0 * k1[coordIndex]
						  + b10_1 * k2[coordIndex]
						  + b10_2 * k3[coordIndex]
						  + b10_3 * k4[coordIndex]
						  + b10_4 * k5[coordIndex]
						  + b10_5 * k6[coordIndex]
						  + b10_6 * k7[coordIndex]
						  + b10_7 * k8[coordIndex]
						  + b10_8 * k9[coordIndex]
						  + b10_9 * k10[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep *a9, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 10:
		{
				   Scalar a10(0.700701039770150737151099854830749337941407049265546408969222);
				   Scalar b11_0(0.0419700073362782579861792864787277787213483656543104611245994);
				   Scalar b11_1(0);
				   Scalar b11_2(0);
				   Scalar b11_3(0);
				   Scalar b11_4(0);
				   Scalar b11_5(-0.00317987696266205093901912847692712407988609169703103952205634);
				   Scalar b11_6(0.806397714906192077260821711520379506393543111567419750119748);
				   Scalar b11_7(0.0975983126412388979093522850684288851314672048003054550357187);
				   Scalar b11_8(0.778575578158398909027512446452927238999763460594181964958853);
				   Scalar b11_9(0.204890423831599428189499202098105603312029235081420653574829);
				   Scalar b11_10(-1.56261579627468188307070943950527825211462892236424360892806);
				   system->GetCurrDerivatives(k11, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b11_0 * k1[coordIndex]
						   + b11_1 * k2[coordIndex]
						   + b11_2 * k3[coordIndex]
						   + b11_3 * k4[coordIndex]
						   + b11_4 * k5[coordIndex]
						   + b11_5 * k6[coordIndex]
						   + b11_6 * k7[coordIndex]
						   + b11_7 * k8[coordIndex]
						   + b11_8 * k9[coordIndex]
						   + b11_9 * k10[coordIndex]
						   + b11_10 * k11[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a10, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 11:
		{
				   Scalar a11(0.363636363636363636363636363636363636363636363636363636363636);
				   Scalar b12_0(0.0437726782233730163574465242495339811688214967071614123256973);
				   Scalar b12_1(0);
				   Scalar b12_2(0);
				   Scalar b12_3(0);
				   Scalar b12_4(0);
				   Scalar b12_5(0);
				   Scalar b12_6(0);
				   Scalar b12_7(0);
				   Scalar b12_8(0.00624365027520195208794358628580933625281631216903095917201250);
				   Scalar b12_9(0.200043097109577314994435165469647856829066232218264969608768);
				   Scalar b12_10(-0.00805328367804983036823857162048902911923392887337029314844206);
				   Scalar b12_11(0.0211517528067396521915711903523399601316877825157550573051221);
				   system->GetCurrDerivatives(k12, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b12_0 * k1[coordIndex]
						   + b12_1 * k2[coordIndex]
						   + b12_2 * k3[coordIndex]
						   + b12_3 * k4[coordIndex]
						   + b12_4 * k5[coordIndex]
						   + b12_5 * k6[coordIndex]
						   + b12_6 * k7[coordIndex]
						   + b12_7 * k8[coordIndex]
						   + b12_8 * k9[coordIndex]
						   + b12_9 * k10[coordIndex]
						   + b12_10 * k11[coordIndex]
						   + b12_11 * k12[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a11, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 12:
		{
				   Scalar a12(0.263157894736842105263157894736842105263157894736842105263158);
				   Scalar b13_0(0.0283499250363514563095023591920717312247137654896477097768495);
				   Scalar b13_1(0);
				   Scalar b13_2(0);
				   Scalar b13_3(0);
				   Scalar b13_4(0);
				   Scalar b13_5(0);
				   Scalar b13_6(0);
				   Scalar b13_7(0);
				   Scalar b13_8(0.00249163204855817407538949148805995149459884653585417680098222);
				   Scalar b13_9(0.0230138787854593149638399846373742768772087122638142234223658);
				   Scalar b13_10(-0.00322155956692977098724476092467120878189463604760620461043308);
				   Scalar b13_11(0.00988442549447664668946335414487885256040819982786014648129297);
				   Scalar b13_12(-0.0213010771328887351384307642875927384886634565429572466632092);
				   system->GetCurrDerivatives(k13, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b13_0 * k1[coordIndex]
						   + b13_1 * k2[coordIndex]
						   + b13_2 * k3[coordIndex]
						   + b13_3 * k4[coordIndex]
						   + b13_4 * k5[coordIndex]
						   + b13_5 * k6[coordIndex]
						   + b13_6 * k7[coordIndex]
						   + b13_7 * k8[coordIndex]
						   + b13_8 * k9[coordIndex]
						   + b13_9 * k10[coordIndex]
						   + b13_10 * k11[coordIndex]
						   + b13_11 * k12[coordIndex]
						   + b13_12 * k13[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a12, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 13:
		{
				   Scalar a13(0.0392172246650270859125196642501208648863714315266128052078483);
				   Scalar b14_0(0.343511894290243001049432234735147943083353174980701426268122);
				   Scalar b14_1(0);
				   Scalar b14_2(0);
				   Scalar b14_3(0);
				   Scalar b14_4(0);
				   Scalar b14_5(0);
				   Scalar b14_6(0);
				   Scalar b14_7(0);
				   Scalar b14_8(0.210451912023627385609097011999010655788807405225626700040882);
				   Scalar b14_9(1.03427452057230411936482926828825709938667999698324740166559);
				   Scalar b14_10(0.00600303645864422487051240448206640574939078092406156945568306);
				   Scalar b14_11(0.855938125099619537578012106002407728915062652616416005816477);
				   Scalar b14_12(-0.977235005036766810872264852372525633013107656892839677696022);
				   Scalar b14_13(-0.660026980479294694616225013856327693720573981219974874776419);
				   system->GetCurrDerivatives(k14, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b14_0 * k1[coordIndex]
						   + b14_1 * k2[coordIndex]
						   + b14_2 * k3[coordIndex]
						   + b14_3 * k4[coordIndex]
						   + b14_4 * k5[coordIndex]
						   + b14_5 * k6[coordIndex]
						   + b14_6 * k7[coordIndex]
						   + b14_7 * k8[coordIndex]
						   + b14_8 * k9[coordIndex]
						   + b14_9 * k10[coordIndex]
						   + b14_10 * k11[coordIndex]
						   + b14_11 * k12[coordIndex]
						   + b14_12 * k13[coordIndex]
						   + b14_13 * k14[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a13, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 14:
		{
				   Scalar a14(0.812917502928376762983393159278036506189612372617238550774312);
				   Scalar b15_0(-0.0143574001672168069538206399935076366657755954378399880691949);
				   Scalar b15_1(0);
				   Scalar b15_2(0);
				   Scalar b15_3(0);
				   Scalar b15_4(0);
				   Scalar b15_5(0);
				   Scalar b15_6(0);
				   Scalar b15_7(0);
				   Scalar b15_8(-0.0366253270049039970293685796848974791733119081733552207318285);
				   Scalar b15_9(0.0350254975636213681976849406979846524346789082471103574920148);
				   Scalar b15_10(0.0360946016362113508931786658758335239823689929864237671348749);
				   Scalar b15_11(-0.0265219967553681106351595946834601923649627012457464284442911);
				   Scalar b15_12(0.0445699011305698119638911537508839908104336323082226770910408);
				   Scalar b15_13(0.124343093331358243286225595741786448038973408895106741855721);
				   Scalar b15_14(0.00413829693239480694403512496204335960426192908674476033832967);
				   system->GetCurrDerivatives(k15, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b15_0 * k1[coordIndex]
						   + b15_1 * k2[coordIndex]
						   + b15_2 * k3[coordIndex]
						   + b15_3 * k4[coordIndex]
						   + b15_4 * k5[coordIndex]
						   + b15_5 * k6[coordIndex]
						   + b15_6 * k7[coordIndex]
						   + b15_7 * k8[coordIndex]
						   + b15_8 * k9[coordIndex]
						   + b15_9 * k10[coordIndex]
						   + b15_10 * k11[coordIndex]
						   + b15_11 * k12[coordIndex]
						   + b15_12 * k13[coordIndex]
						   + b15_13 * k14[coordIndex]
						   + b15_14 * k15[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a14, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 15:
		{
				   Scalar a15(0.166666666666666666666666666666666666666666666666666666666667);
				   Scalar b16_0(0.356032404425120290975609116398089176264106222379748802654822);
				   Scalar b16_1(0);
				   Scalar b16_2(0);
				   Scalar b16_3(0);
				   Scalar b16_4(0);
				   Scalar b16_5(0);
				   Scalar b16_6(0);
				   Scalar b16_7(0);
				   Scalar b16_8(-0.450192758947562595966821779075956175110645100214763601190349);
				   Scalar b16_9(0.430527907083710898626656292808782917793030154094709462877146);
				   Scalar b16_10(0.511973029011022237668556960394071692077125787030651386389972);
				   Scalar b16_11(0.908303638886404260390159124638110213997496214819904630546596);
				   Scalar b16_12(-1.23921093371933931757372469151534028854413889248605726186520);
				   Scalar b16_13(-0.649048661671761465141672348879062553905402831967191097656668);
				   Scalar b16_14(0.251708904586819292210480529948970541404887852931447491219418);
				   Scalar b16_15(0.779906470345586398810756795282334476023540593411550187024263);
				   system->GetCurrDerivatives(k16, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b16_0 * k1[coordIndex]
						   + b16_1 * k2[coordIndex]
						   + b16_2 * k3[coordIndex]
						   + b16_3 * k4[coordIndex]
						   + b16_4 * k5[coordIndex]
						   + b16_5 * k6[coordIndex]
						   + b16_6 * k7[coordIndex]
						   + b16_7 * k8[coordIndex]
						   + b16_8 * k9[coordIndex]
						   + b16_9 * k10[coordIndex]
						   + b16_10 * k11[coordIndex]
						   + b16_11 * k12[coordIndex]
						   + b16_12 * k13[coordIndex]
						   + b16_13 * k14[coordIndex]
						   + b16_14 * k15[coordIndex]
						   + b16_15 * k16[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a15, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 16:
		{
				   Scalar a16(0.900000000000000000000000000000000000000000000000000000000000);
				   Scalar b17_0(0.0130935687406513066406881206418834980127470438213192487844956);
				   Scalar b17_1(0);
				   Scalar b17_2(0);
				   Scalar b17_3(0);
				   Scalar b17_4(0);
				   Scalar b17_5(0);
				   Scalar b17_6(0);
				   Scalar b17_7(0);
				   Scalar b17_8(0);
				   Scalar b17_9(0);
				   Scalar b17_10(0);
				   Scalar b17_11(0);
				   Scalar b17_12(-0.0000932053067985113945908461962767108237858631509684667142124826);
				   Scalar b17_13(0.0505374334262299359640090443138590726770942344716122381702746);
				   Scalar b17_14(8.04470341944487979109579109610197797641311868930865361048975E-7);
				   Scalar b17_15(0.000591726029494171190528755742777717259844340971924321528178248);
				   Scalar b17_16(-4.01614722154557337064691684906375587732264247950093804676867E-7);
				   system->GetCurrDerivatives(k17, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b17_0 * k1[coordIndex]
						   + b17_1 * k2[coordIndex]
						   + b17_2 * k3[coordIndex]
						   + b17_3 * k4[coordIndex]
						   + b17_4 * k5[coordIndex]
						   + b17_5 * k6[coordIndex]
						   + b17_6 * k7[coordIndex]
						   + b17_7 * k8[coordIndex]
						   + b17_8 * k9[coordIndex]
						   + b17_9 * k10[coordIndex]
						   + b17_10 * k11[coordIndex]
						   + b17_11 * k12[coordIndex]
						   + b17_12 * k13[coordIndex]
						   + b17_13 * k14[coordIndex]
						   + b17_14 * k15[coordIndex]
						   + b17_15 * k16[coordIndex]
						   + b17_16 * k17[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a16, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 17:
		{
				   Scalar a17(0.0641299257451966923312771193896682809481096651615083225402924);
				   Scalar b18_0(0.0207926484466053012541944544000765652167255206144373407979758);
				   Scalar b18_1(0);
				   Scalar b18_2(0);
				   Scalar b18_3(0);
				   Scalar b18_4(0);
				   Scalar b18_5(0);
				   Scalar b18_6(0);
				   Scalar b18_7(0);
				   Scalar b18_8(0);
				   Scalar b18_9(0);
				   Scalar b18_10(0);
				   Scalar b18_11(0);
				   Scalar b18_12(0.000582695918800085915101902697837284108951406103029871570103075);
				   Scalar b18_13(-0.00801700732358815939083342186525852746640558465919633524655451);
				   Scalar b18_14(4.03847643847136940375170821743560570484117290330895506618968E-6);
				   Scalar b18_15(0.0854609998055506144225056114567535602510114622033622491802597);
				   Scalar b18_16(-2.04486480935804242706707569691004307904442837552677456232848E-6);
				   Scalar b18_17(0.105328578824431893399799402979093997354240904235172843146582);
				   system->GetCurrDerivatives(k18, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b18_0 * k1[coordIndex]
						   + b18_1 * k2[coordIndex]
						   + b18_2 * k3[coordIndex]
						   + b18_3 * k4[coordIndex]
						   + b18_4 * k5[coordIndex]
						   + b18_5 * k6[coordIndex]
						   + b18_6 * k7[coordIndex]
						   + b18_7 * k8[coordIndex]
						   + b18_8 * k9[coordIndex]
						   + b18_9 * k10[coordIndex]
						   + b18_10 * k11[coordIndex]
						   + b18_11 * k12[coordIndex]
						   + b18_12 * k13[coordIndex]
						   + b18_13 * k14[coordIndex]
						   + b18_14 * k15[coordIndex]
						   + b18_15 * k16[coordIndex]
						   + b18_16 * k17[coordIndex]
						   + b18_17 * k18[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a17, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 18:
		{
				   Scalar a18(0.204149909283428848927744634301023405027149505241333751628870);
				   Scalar b19_0(1.40153449795736021415446247355771306718486452917597731683689);
				   Scalar b19_1(0);
				   Scalar b19_2(0);
				   Scalar b19_3(0);
				   Scalar b19_4(0);
				   Scalar b19_5(0);
				   Scalar b19_6(0);
				   Scalar b19_7(0);
				   Scalar b19_8(0);
				   Scalar b19_9(0);
				   Scalar b19_10(0);
				   Scalar b19_11(0);
				   Scalar b19_12(-0.230252000984221261616272410367415621261130298274455611733277);
				   Scalar b19_13(-7.21106840466912905659582237106874247165856493509961561958267);
				   Scalar b19_14(0.00372901560694836335236995327852132340217759566678662385552634);
				   Scalar b19_15(-4.71415495727125020678778179392224757011323373221820091641216);
				   Scalar b19_16(-0.00176367657545349242053841995032797673574903886695600132759652);
				   Scalar b19_17(7.64130548038698765563029310880237651185173367813936997648198);
				   Scalar b19_18(3.50602043659751834989896082949744710968212949893375368243588);
				   system->GetCurrDerivatives(k19, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b19_0 * k1[coordIndex]
						   + b19_1 * k2[coordIndex]
						   + b19_2 * k3[coordIndex]
						   + b19_3 * k4[coordIndex]
						   + b19_4 * k5[coordIndex]
						   + b19_5 * k6[coordIndex]
						   + b19_6 * k7[coordIndex]
						   + b19_7 * k8[coordIndex]
						   + b19_8 * k9[coordIndex]
						   + b19_9 * k10[coordIndex]
						   + b19_10 * k11[coordIndex]
						   + b19_11 * k12[coordIndex]
						   + b19_12 * k13[coordIndex]
						   + b19_13 * k14[coordIndex]
						   + b19_14 * k15[coordIndex]
						   + b19_15 * k16[coordIndex]
						   + b19_16 * k17[coordIndex]
						   + b19_17 * k18[coordIndex]
						   + b19_18 * k19[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a18, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 19:
		{
				   Scalar a19(0.395350391048760565615671369827324372352227297456659450554577);
				   Scalar b20_0(11.9514650694120686799372385830716401674473610826553517297976);
				   Scalar b20_1(0);
				   Scalar b20_2(0);
				   Scalar b20_3(0);
				   Scalar b20_4(0);
				   Scalar b20_5(0);
				   Scalar b20_6(0);
				   Scalar b20_7(0);
				   Scalar b20_8(0);
				   Scalar b20_9(0);
				   Scalar b20_10(0);
				   Scalar b20_11(0);
				   Scalar b20_12(7.79480932108175968783516700231764388220284279598980948538579);
				   Scalar b20_13(-56.4501393867325792523560991120904281440468100061340556540132);
				   Scalar b20_14(0.0912376306930644901344530449290276645709607450403673704844997);
				   Scalar b20_15(-12.7336279925434886201945524309199275038162717529918963305155);
				   Scalar b20_16(-0.0396895921904719712313542810939736674712383070433147873009352);
				   Scalar b20_17(54.4392141883570886996225765155307791861438378423305337073797);
				   Scalar b20_18(-3.64411637921569236846406990361350645806721478409266709351203);
				   Scalar b20_19(-0.804503249910509910899030787958579499315694913210787878260459);
				   system->GetCurrDerivatives(k20, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b20_0 * k1[coordIndex]
						   + b20_1 * k2[coordIndex]
						   + b20_2 * k3[coordIndex]
						   + b20_3 * k4[coordIndex]
						   + b20_4 * k5[coordIndex]
						   + b20_5 * k6[coordIndex]
						   + b20_6 * k7[coordIndex]
						   + b20_7 * k8[coordIndex]
						   + b20_8 * k9[coordIndex]
						   + b20_9 * k10[coordIndex]
						   + b20_10 * k11[coordIndex]
						   + b20_11 * k12[coordIndex]
						   + b20_12 * k13[coordIndex]
						   + b20_13 * k14[coordIndex]
						   + b20_14 * k15[coordIndex]
						   + b20_15 * k16[coordIndex]
						   + b20_16 * k17[coordIndex]
						   + b20_17 * k18[coordIndex]
						   + b20_18 * k19[coordIndex]
						   + b20_19 * k20[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a19, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 20:
		{
				   Scalar a20(0.604649608951239434384328630172675627647772702543340549445423);
				   Scalar b21_0(-148.809426507100488427838868268647625561930612082148597076690);
				   Scalar b21_1(0);
				   Scalar b21_2(0);
				   Scalar b21_3(0);
				   Scalar b21_4(0);
				   Scalar b21_5(0);
				   Scalar b21_6(0);
				   Scalar b21_7(0);
				   Scalar b21_8(0);
				   Scalar b21_9(0);
				   Scalar b21_10(0);
				   Scalar b21_11(0);
				   Scalar b21_12(-91.7295278291256484357935662402321623495228729036354276506427);
				   Scalar b21_13(707.656144971598359834575719286335716154821128966649565194286);
				   Scalar b21_14(-1.10563611857482440905296961311590930801338308942637769555540);
				   Scalar b21_15(176.134591883811372587859898076055660406999516762301689616841);
				   Scalar b21_16(0.491384824214880662268898345164454557416884631402764792538746);
				   Scalar b21_17(-684.278000449814944358237535610895081956077167893600278300805);
				   Scalar b21_18(27.9910604998398258984224332124380407446002518400668657974589);
				   Scalar b21_19(13.1939710030282333443670964371153238435064159623744975073252);
				   Scalar b21_20(1.25128781283980445450114974148056006317268830077396406361417);
				   system->GetCurrDerivatives(k21, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b21_0 * k1[coordIndex]
						   + b21_1 * k2[coordIndex]
						   + b21_2 * k3[coordIndex]
						   + b21_3 * k4[coordIndex]
						   + b21_4 * k5[coordIndex]
						   + b21_5 * k6[coordIndex]
						   + b21_6 * k7[coordIndex]
						   + b21_7 * k8[coordIndex]
						   + b21_8 * k9[coordIndex]
						   + b21_9 * k10[coordIndex]
						   + b21_10 * k11[coordIndex]
						   + b21_11 * k12[coordIndex]
						   + b21_12 * k13[coordIndex]
						   + b21_13 * k14[coordIndex]
						   + b21_14 * k15[coordIndex]
						   + b21_15 * k16[coordIndex]
						   + b21_16 * k17[coordIndex]
						   + b21_17 * k18[coordIndex]
						   + b21_18 * k19[coordIndex]
						   + b21_19 * k20[coordIndex]
						   + b21_20 * k21[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a20, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 21:
		{
				   Scalar a21(0.795850090716571151072255365698976594972850494758666248371130);
				   Scalar b22_0(-9.67307946948196763644126118433219395839951408571877262880482);
				   Scalar b22_1(0);
				   Scalar b22_2(0);
				   Scalar b22_3(0);
				   Scalar b22_4(0);
				   Scalar b22_5(0);
				   Scalar b22_6(0);
				   Scalar b22_7(0);
				   Scalar b22_8(0);
				   Scalar b22_9(0);
				   Scalar b22_10(0);
				   Scalar b22_11(0);
				   Scalar b22_12(-4.46990150858505531443846227701960360497830681408751431146712);
				   Scalar b22_13(45.5127128690952681968241950400052751178905907817398483534845);
				   Scalar b22_14(-0.0713085086183826912791492024438246129930559805352394367050813);
				   Scalar b22_15(11.2273614068412741582590624479939384207826800776794485051540);
				   Scalar b22_16(0.126244376717622724516237912909138809361786889819105426371393);
				   Scalar b22_17(-43.5439339549483313605810624907242107623814304467621407753424);
				   Scalar b22_18(0.787174307543058978398792994996550902064546091443233850464377);
				   Scalar b22_19(0.532264696744684215669300708603886690785395776821503851830821);
				   Scalar b22_20(0.422422733996325326010225127471388772575086538809603346825334);
				   Scalar b22_21(0.0859131249503067107308438031499859443441115056294154956487671);
				   system->GetCurrDerivatives(k22, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b22_0 * k1[coordIndex]
						   + b22_1 * k2[coordIndex]
						   + b22_2 * k3[coordIndex]
						   + b22_3 * k4[coordIndex]
						   + b22_4 * k5[coordIndex]
						   + b22_5 * k6[coordIndex]
						   + b22_6 * k7[coordIndex]
						   + b22_7 * k8[coordIndex]
						   + b22_8 * k9[coordIndex]
						   + b22_9 * k10[coordIndex]
						   + b22_10 * k11[coordIndex]
						   + b22_11 * k12[coordIndex]
						   + b22_12 * k13[coordIndex]
						   + b22_13 * k14[coordIndex]
						   + b22_14 * k15[coordIndex]
						   + b22_15 * k16[coordIndex]
						   + b22_16 * k17[coordIndex]
						   + b22_17 * k18[coordIndex]
						   + b22_18 * k19[coordIndex]
						   + b22_19 * k20[coordIndex]
						   + b22_20 * k21[coordIndex]
						   + b22_21 * k22[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a21, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 22:
		{
				   Scalar a22(0.935870074254803307668722880610331719051890334838491677459708);
				   Scalar b23_0(-10.0664032447054702403396606900426891472202824757968765569183);
				   Scalar b23_1(0);
				   Scalar b23_2(0);
				   Scalar b23_3(0);
				   Scalar b23_4(0);
				   Scalar b23_5(0);
				   Scalar b23_6(0);
				   Scalar b23_7(0);
				   Scalar b23_8(-0.0366253270049039970293685796848974791733119081733552207318285);
				   Scalar b23_9(0.0350254975636213681976849406979846524346789082471103574920148);
				   Scalar b23_10(0.0360946016362113508931786658758335239823689929864237671348749);
				   Scalar b23_11(-0.0265219967553681106351595946834601923649627012457464284442911);
				   Scalar b23_12(-6.27088972181464143590553149478871603839356122957396018530209);
				   Scalar b23_13(48.2079237442562989090702103008195063923492593141636117832993);
				   Scalar b23_14(-0.0694471689136165640882395180583732834557754169149088630301342);
				   Scalar b23_15(12.6810690204850295698341370913609807066108483811412127009785);
				   Scalar b23_16(0.0119671168968323754838161435501011294100927813964199613229864);
				   Scalar b23_17(-46.7249764992482408003358268242662695593201321659795608950429);
				   Scalar b23_18(1.33029613326626711314710039298216591399033511191227101321435);
				   Scalar b23_19(1.00766787503398298353438903619926657771162717793661719708370);
				   Scalar b23_20(0.0209512051933665091664122388475480702892770753864487241177616);
				   Scalar b23_21(0.0210134706331264177317735424331396407424412188443757490871603);
				   Scalar b23_22(0.00952196014417121794175101542454575907376360233658356240547761);
				   system->GetCurrDerivatives(k23, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b23_0 * k1[coordIndex]
						   + b23_1 * k2[coordIndex]
						   + b23_2 * k3[coordIndex]
						   + b23_3 * k4[coordIndex]
						   + b23_4 * k5[coordIndex]
						   + b23_5 * k6[coordIndex]
						   + b23_6 * k7[coordIndex]
						   + b23_7 * k8[coordIndex]
						   + b23_8 * k9[coordIndex]
						   + b23_9 * k10[coordIndex]
						   + b23_10 * k11[coordIndex]
						   + b23_11 * k12[coordIndex]
						   + b23_12 * k13[coordIndex]
						   + b23_13 * k14[coordIndex]
						   + b23_14 * k15[coordIndex]
						   + b23_15 * k16[coordIndex]
						   + b23_16 * k17[coordIndex]
						   + b23_17 * k18[coordIndex]
						   + b23_18 * k19[coordIndex]
						   + b23_19 * k20[coordIndex]
						   + b23_20 * k21[coordIndex]
						   + b23_21 * k22[coordIndex]
						   + b23_22 * k23[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a22, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 23:
		{
				   Scalar a23(0.166666666666666666666666666666666666666666666666666666666667);
				   Scalar b24_0(-409.478081677743708772589097409370357624424341606752069725341);
				   Scalar b24_1(0);
				   Scalar b24_2(0);
				   Scalar b24_3(0);
				   Scalar b24_4(0);
				   Scalar b24_5(0);
				   Scalar b24_6(0);
				   Scalar b24_7(0);
				   Scalar b24_8(0.210451912023627385609097011999010655788807405225626700040882);
				   Scalar b24_9(1.03427452057230411936482926828825709938667999698324740166559);
				   Scalar b24_10(0.00600303645864422487051240448206640574939078092406156945568306);
				   Scalar b24_11(0.855938125099619537578012106002407728915062652616416005816477);
				   Scalar b24_12(-250.516998547447860492777657729316130386584050420782075966990);
				   Scalar b24_13(1946.42466652388427766053750328264758595829850895761428240231);
				   Scalar b24_14(-3.04503882102310365506105809086860882786950544097602101685174);
				   Scalar b24_15(490.626379528281713521208265299168083841598542274061671576230);
				   Scalar b24_16(1.56647589531270907115484067013597445739595615245966775329993);
				   Scalar b24_17(-1881.97428994011173362217267377035870619215906638453056643641);
				   Scalar b24_18(75.2592224724847175278837713643303149821620618914245864351135);
				   Scalar b24_19(34.5734356980331067622434344736554689696728644793551014989002);
				   Scalar b24_20(3.21147679440968961435417361847073755169022966748891627882572);
				   Scalar b24_21(-0.460408041738414391307201404237058848867245095265382820823055);
				   Scalar b24_22(-0.0870718339841810522431884137957986245724252047388936572215438);
				   Scalar b24_23(-7.39351814158303067567016952195521063999185773249132944724553);
				   system->GetCurrDerivatives(k24, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b24_0 * k1[coordIndex]
						   + b24_1 * k2[coordIndex]
						   + b24_2 * k3[coordIndex]
						   + b24_3 * k4[coordIndex]
						   + b24_4 * k5[coordIndex]
						   + b24_5 * k6[coordIndex]
						   + b24_6 * k7[coordIndex]
						   + b24_7 * k8[coordIndex]
						   + b24_8 * k9[coordIndex]
						   + b24_9 * k10[coordIndex]
						   + b24_10 * k11[coordIndex]
						   + b24_11 * k12[coordIndex]
						   + b24_12 * k13[coordIndex]
						   + b24_13 * k14[coordIndex]
						   + b24_14 * k15[coordIndex]
						   + b24_15 * k16[coordIndex]
						   + b24_16 * k17[coordIndex]
						   + b24_17 * k18[coordIndex]
						   + b24_18 * k19[coordIndex]
						   + b24_19 * k20[coordIndex]
						   + b24_20 * k21[coordIndex]
						   + b24_21 * k22[coordIndex]
						   + b24_22 * k23[coordIndex]
						   + b24_23 * k24[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a23, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 24:
		{
				   Scalar a24(0.812917502928376762983393159278036506189612372617238550774312);
				   Scalar b25_0(3.43347475853550878921093496257596781120623891072008459930197);
				   Scalar b25_1(0);
				   Scalar b25_2(0);
				   Scalar b25_3(0);
				   Scalar b25_4(0);
				   Scalar b25_5(0);
				   Scalar b25_6(0);
				   Scalar b25_7(0);
				   Scalar b25_8(0.00249163204855817407538949148805995149459884653585417680098222);
				   Scalar b25_9(0.0230138787854593149638399846373742768772087122638142234223658);
				   Scalar b25_10(-0.00322155956692977098724476092467120878189463604760620461043308);
				   Scalar b25_11(0.00988442549447664668946335414487885256040819982786014648129297);
				   Scalar b25_12(2.16252799377922507788307841904757354045759225335732707916530);
				   Scalar b25_13(-16.2699864546457421328065640660139489006987552040228852402716);
				   Scalar b25_14(-0.128534502120524552843583417470935010538029037542654506231743);
				   Scalar b25_15(-8.98915042666504253089307820833379330486511746063552853023189);
				   Scalar b25_16(-0.00348595363232025333387080201851013650192401767250513765000963);
				   Scalar b25_17(15.7936194113339807536235187388695574135853387025139738341334);
				   Scalar b25_18(-0.574403330914095065628165482017335820148383663195675408024658);
				   Scalar b25_19(-0.345602039021393296692722496608124982535237228827655306030152);
				   Scalar b25_20(-0.00662241490206585091731619991383757781133067992707418687587487);
				   Scalar b25_21(-0.00777788129242204164032546458607364309759347209626759111946150);
				   Scalar b25_22(-0.00356084192402274913338827232697437364675240818791706587952939);
				   Scalar b25_23(4.79282506449930799649797749629840189457296934139359048988332);
				   Scalar b25_24(0.153725464873068577844576387402512082757034273069877432944621);
				   system->GetCurrDerivatives(k25, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b25_0 * k1[coordIndex]
						   + b25_1 * k2[coordIndex]
						   + b25_2 * k3[coordIndex]
						   + b25_3 * k4[coordIndex]
						   + b25_4 * k5[coordIndex]
						   + b25_5 * k6[coordIndex]
						   + b25_6 * k7[coordIndex]
						   + b25_7 * k8[coordIndex]
						   + b25_8 * k9[coordIndex]
						   + b25_9 * k10[coordIndex]
						   + b25_10 * k11[coordIndex]
						   + b25_11 * k12[coordIndex]
						   + b25_12 * k13[coordIndex]
						   + b25_13 * k14[coordIndex]
						   + b25_14 * k15[coordIndex]
						   + b25_15 * k16[coordIndex]
						   + b25_16 * k17[coordIndex]
						   + b25_17 * k18[coordIndex]
						   + b25_18 * k19[coordIndex]
						   + b25_19 * k20[coordIndex]
						   + b25_20 * k21[coordIndex]
						   + b25_21 * k22[coordIndex]
						   + b25_22 * k23[coordIndex]
						   + b25_23 * k24[coordIndex]
						   + b25_24 * k25[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a24, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 25:
		{
				   Scalar a25(0.0392172246650270859125196642501208648863714315266128052078483);
				   Scalar b26_0(32.3038520871985442326994734440031535091364975047784630088983);
				   Scalar b26_1(0);
				   Scalar b26_2(0);
				   Scalar b26_3(0);
				   Scalar b26_4(0);
				   Scalar b26_5(-0.00317987696266205093901912847692712407988609169703103952205634);
				   Scalar b26_6(0.806397714906192077260821711520379506393543111567419750119748);
				   Scalar b26_7(0.0975983126412388979093522850684288851314672048003054550357187);
				   Scalar b26_8(0.778575578158398909027512446452927238999763460594181964958853);
				   Scalar b26_9(0.204890423831599428189499202098105603312029235081420653574829);
				   Scalar b26_10(-1.56261579627468188307070943950527825211462892236424360892806);
				   Scalar b26_11(0);
				   Scalar b26_12(16.3429891882310570648504243973927174708753353504154550405647);
				   Scalar b26_13(-154.544555293543621230730189631471036399316683669609116705323);
				   Scalar b26_14(1.56971088703334872692034283417621761466263593582497085955201);
				   Scalar b26_15(3.27685545087248131321429817269900731165522404974733504794135);
				   Scalar b26_16(-0.0503489245193653176348040727199783626534081095691632396802451);
				   Scalar b26_17(153.321151858041665070593767885914694011224363102594556731397);
				   Scalar b26_18(7.17568186327720495846766484814784143567826308034865369443637);
				   Scalar b26_19(-2.94036748675300481945917659896930989215320594380777597403592);
				   Scalar b26_20(-0.0665845946076803144470749676022628870281920493197256887985612);
				   Scalar b26_21(-0.0462346054990843661229248668562217261176966514016859284197145);
				   Scalar b26_22(-0.0204198733585679401539388228617269778848579774821581777675337);
				   Scalar b26_23(-53.3523106438735850515953441165998107974045090495791591218714);
				   Scalar b26_24(-1.35548714715078654978732186705996404017554501614191325114947);
				   Scalar b26_25(-1.57196275801232751882901735171459249177687219114442583461866);
				   system->GetCurrDerivatives(k26, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b26_0 * k1[coordIndex]
						   + b26_1 * k2[coordIndex]
						   + b26_2 * k3[coordIndex]
						   + b26_3 * k4[coordIndex]
						   + b26_4 * k5[coordIndex]
						   + b26_5 * k6[coordIndex]
						   + b26_6 * k7[coordIndex]
						   + b26_7 * k8[coordIndex]
						   + b26_8 * k9[coordIndex]
						   + b26_9 * k10[coordIndex]
						   + b26_10 * k11[coordIndex]
						   + b26_11 * k12[coordIndex]
						   + b26_12 * k13[coordIndex]
						   + b26_13 * k14[coordIndex]
						   + b26_14 * k15[coordIndex]
						   + b26_15 * k16[coordIndex]
						   + b26_16 * k17[coordIndex]
						   + b26_17 * k18[coordIndex]
						   + b26_18 * k19[coordIndex]
						   + b26_19 * k20[coordIndex]
						   + b26_20 * k21[coordIndex]
						   + b26_21 * k22[coordIndex]
						   + b26_22 * k23[coordIndex]
						   + b26_23 * k24[coordIndex]
						   + b26_24 * k25[coordIndex]
						   + b26_25 * k26[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a25, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 26:
		{
				   Scalar a26(0.363636363636363636363636363636363636363636363636363636363636);
				   Scalar b27_0(-16.6451467486341512872031294403931758764560371130818978459405);
				   Scalar b27_1(0);
				   Scalar b27_2(0);
				   Scalar b27_3(0);
				   Scalar b27_4(0);
				   Scalar b27_5(0.00592232780324503308042990005798046524738389560444257136834990);
				   Scalar b27_6(0.470326159963841112217224303205894113455362530746108825010848);
				   Scalar b27_7(0.299688863848679000853981837096192399136831121671781279184194);
				   Scalar b27_8(-0.247656877593994914689992276329810825853958069263947095548189);
				   Scalar b27_9(0.110895029771437682893999851839061714522445173600678718208625);
				   Scalar b27_10(0);
				   Scalar b27_11(-0.491719043846229147070666628704194097678081907210673044988866);
				   Scalar b27_12(-11.4743154427289496968389492564352536350842454130853175250727);
				   Scalar b27_13(80.2593166576230272541702485886484400152793366623589989106256);
				   Scalar b27_14(-0.384132303980042847625312526759029103746926841342088219165648);
				   Scalar b27_15(7.28147667468107583471326950926136115767612581862877764249646);
				   Scalar b27_16(-0.132699384612248379510571708176035274836827341616751884314074);
				   Scalar b27_17(-81.0799832525730726674679289752255240006070716633632990308935);
				   Scalar b27_18(-1.25037492835620639521768185656179119962253747492403205797494);
				   Scalar b27_19(2.59263594969543681023776379504377324994226447359296887778718);
				   Scalar b27_20(-0.301440298346404539830163997260526875264431537275641495291993);
				   Scalar b27_21(0.221384460789832337451706451572773791695246839057318414301020);
				   Scalar b27_22(0.0827577274771892931955989870974693152996276435429809890551210);
				   Scalar b27_23(18.9960662040611520464672450037243263998175161412237156872211);
				   Scalar b27_24(0.269231946409639685623468015128334167460051910348912845121977);
				   Scalar b27_25(1.62674827447066537462989364929628933988125029284183680279020);
				   Scalar b27_26(0.491719043846229147070666628704194097678081907210673044988866);
				   system->GetCurrDerivatives(k27, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b27_0 * k1[coordIndex]
						   + b27_1 * k2[coordIndex]
						   + b27_2 * k3[coordIndex]
						   + b27_3 * k4[coordIndex]
						   + b27_4 * k5[coordIndex]
						   + b27_5 * k6[coordIndex]
						   + b27_6 * k7[coordIndex]
						   + b27_7 * k8[coordIndex]
						   + b27_8 * k9[coordIndex]
						   + b27_9 * k10[coordIndex]
						   + b27_10 * k11[coordIndex]
						   + b27_11 * k12[coordIndex]
						   + b27_12 * k13[coordIndex]
						   + b27_13 * k14[coordIndex]
						   + b27_14 * k15[coordIndex]
						   + b27_15 * k16[coordIndex]
						   + b27_16 * k17[coordIndex]
						   + b27_17 * k18[coordIndex]
						   + b27_18 * k19[coordIndex]
						   + b27_19 * k20[coordIndex]
						   + b27_20 * k21[coordIndex]
						   + b27_21 * k22[coordIndex]
						   + b27_22 * k23[coordIndex]
						   + b27_23 * k24[coordIndex]
						   + b27_24 * k25[coordIndex]
						   + b27_25 * k26[coordIndex]
						   + b27_26 * k27[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a26, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 27:
		{
				   Scalar a27(0.700701039770150737151099854830749337941407049265546408969222);
				   Scalar b28_0(0.0838479812409052664616968791372814085980533139224911131069335);
				   Scalar b28_1(0);
				   Scalar b28_2(0);
				   Scalar b28_3(0);
				   Scalar b28_4(0);
				   Scalar b28_5(-0.0117949367100973814319755056031295775367961960590736150777613);
				   Scalar b28_6(-0.247299020568812652339473838743194598325992840353340132697498);
				   Scalar b28_7(0.0978080858367729012259313014081291665503740655476733940756599);
				   Scalar b28_8(0.217590689243420631360008651767860318344168120024782176879989);
				   Scalar b28_9(0);
				   Scalar b28_10(0.137585606763325224865659632196787746647447222975084865975440);
				   Scalar b28_11(0.0439870229715046685058790092341545026046103890294261359042581);
				   Scalar b28_12(0);
				   Scalar b28_13(-0.513700813768193341957004456618630303738757363641964030086972);
				   Scalar b28_14(0.826355691151315508644211308399153458701423158616168576922372);
				   Scalar b28_15(25.7018139719811832625873882972519939511136556341960074626615);
				   Scalar b28_16(0);
				   Scalar b28_17(0);
				   Scalar b28_18(0);
				   Scalar b28_19(0);
				   Scalar b28_20(0);
				   Scalar b28_21(0);
				   Scalar b28_22(0);
				   Scalar b28_23(-25.7018139719811832625873882972519939511136556341960074626615);
				   Scalar b28_24(-0.826355691151315508644211308399153458701423158616168576922372);
				   Scalar b28_25(0.513700813768193341957004456618630303738757363641964030086972);
				   Scalar b28_26(-0.0439870229715046685058790092341545026046103890294261359042581);
				   Scalar b28_27(-0.137585606763325224865659632196787746647447222975084865975440);
				   system->GetCurrDerivatives(k28, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b28_0 * k1[coordIndex]
						   + b28_1 * k2[coordIndex]
						   + b28_2 * k3[coordIndex]
						   + b28_3 * k4[coordIndex]
						   + b28_4 * k5[coordIndex]
						   + b28_5 * k6[coordIndex]
						   + b28_6 * k7[coordIndex]
						   + b28_7 * k8[coordIndex]
						   + b28_8 * k9[coordIndex]
						   + b28_9 * k10[coordIndex]
						   + b28_10 * k11[coordIndex]
						   + b28_11 * k12[coordIndex]
						   + b28_12 * k13[coordIndex]
						   + b28_13 * k14[coordIndex]
						   + b28_14 * k15[coordIndex]
						   + b28_15 * k16[coordIndex]
						   + b28_16 * k17[coordIndex]
						   + b28_17 * k18[coordIndex]
						   + b28_18 * k19[coordIndex]
						   + b28_19 * k20[coordIndex]
						   + b28_20 * k21[coordIndex]
						   + b28_21 * k22[coordIndex]
						   + b28_22 * k23[coordIndex]
						   + b28_23 * k24[coordIndex]
						   + b28_24 * k25[coordIndex]
						   + b28_25 * k26[coordIndex]
						   + b28_26 * k27[coordIndex]
						   + b28_27 * k28[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a27, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 28:
		{
				   Scalar a28(0.140152799042188765276187487966946717629806463082532936287323);
				   Scalar b29_0(0.124380526654094412881516420868799316268491466359671423163289);
				   Scalar b29_1(0);
				   Scalar b29_2(0);
				   Scalar b29_3(0);
				   Scalar b29_4(0.226120282197584301422238662979202901196752320742633143965145);
				   Scalar b29_5(0.0137885887618080880607695837016477814530969417491493385363543);
				   Scalar b29_6(-0.0672210133996684449749399507414305856950086341525382182856200);
				   Scalar b29_7(0);
				   Scalar b29_8(0);
				   Scalar b29_9(-0.856238975085428354755349769879501772112121597411563802855067);
				   Scalar b29_10(-1.96337522866858908928262850028093813988180440518267404553576);
				   Scalar b29_11(-0.232332822724119401237246257308921847250108199230419994978218);
				   Scalar b29_12(0);
				   Scalar b29_13(4.30660719086453349461668936876562947772432562053478092626764);
				   Scalar b29_14(-2.92722963249465482659787911202390446687687394950633612630592);
				   Scalar b29_15(-82.3131666397858944454492334105458707735761966428138676971041);
				   Scalar b29_16(0);
				   Scalar b29_17(0);
				   Scalar b29_18(0);
				   Scalar b29_19(0);
				   Scalar b29_20(0);
				   Scalar b29_21(0);
				   Scalar b29_22(0);
				   Scalar b29_23(82.3131666397858944454492334105458707735761966428138676971041);
				   Scalar b29_24(2.92722963249465482659787911202390446687687394950633612630592);
				   Scalar b29_25(-4.30660719086453349461668936876562947772432562053478092626764);
				   Scalar b29_26(0.232332822724119401237246257308921847250108199230419994978218);
				   Scalar b29_27(1.96337522866858908928262850028093813988180440518267404553576);
				   Scalar b29_28(0.856238975085428354755349769879501772112121597411563802855067);
				   system->GetCurrDerivatives(k29, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b29_0 * k1[coordIndex]
						   + b29_1 * k2[coordIndex]
						   + b29_2 * k3[coordIndex]
						   + b29_3 * k4[coordIndex]
						   + b29_4 * k5[coordIndex]
						   + b29_5 * k6[coordIndex]
						   + b29_6 * k7[coordIndex]
						   + b29_7 * k8[coordIndex]
						   + b29_8 * k9[coordIndex]
						   + b29_9 * k10[coordIndex]
						   + b29_10 * k11[coordIndex]
						   + b29_11 * k12[coordIndex]
						   + b29_12 * k13[coordIndex]
						   + b29_13 * k14[coordIndex]
						   + b29_14 * k15[coordIndex]
						   + b29_15 * k16[coordIndex]
						   + b29_16 * k17[coordIndex]
						   + b29_17 * k18[coordIndex]
						   + b29_18 * k19[coordIndex]
						   + b29_19 * k20[coordIndex]
						   + b29_20 * k21[coordIndex]
						   + b29_21 * k22[coordIndex]
						   + b29_22 * k23[coordIndex]
						   + b29_23 * k24[coordIndex]
						   + b29_24 * k25[coordIndex]
						   + b29_25 * k26[coordIndex]
						   + b29_26 * k27[coordIndex]
						   + b29_27 * k28[coordIndex]
						   + b29_28 * k29[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a28, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 29:
		{
				   Scalar a29(0.297068384213818357389584716808219413223332094698915687379168);
				   Scalar b30_0(0.103484561636679776672993546511910344499744798201971316606663);
				   Scalar b30_1(0);
				   Scalar b30_2(0);
				   Scalar b30_3(0.122068887306407222589644082868962077139592714834162134741275);
				   Scalar b30_4(0.482574490331246622475134780125688112865919023850168049679402);
				   Scalar b30_5(-0.0381409600015606999730886240005620205664113072478411477421970);
				   Scalar b30_6(0);
				   Scalar b30_7(-0.550499525310802324138388507020508177411414311000037561712836);
				   Scalar b30_8(0);
				   Scalar b30_9(-0.711915811585189227887648262043794387578291882406745570495765);
				   Scalar b30_10(-0.584129605671551340432988730158480872095335329645227595707052);
				   Scalar b30_11(0);
				   Scalar b30_12(0);
				   Scalar b30_13(2.11046308125864932128717300046622750300375054278936987850718);
				   Scalar b30_14(-0.0837494736739572135525742023001037992695260175335123517729291);
				   Scalar b30_15(5.10021499072320914075295969043344113107545060862804249161191);
				   Scalar b30_16(0);
				   Scalar b30_17(0);
				   Scalar b30_18(0);
				   Scalar b30_19(0);
				   Scalar b30_20(0);
				   Scalar b30_21(0);
				   Scalar b30_22(0);
				   Scalar b30_23(-5.10021499072320914075295969043344113107545060862804249161191);
				   Scalar b30_24(0.0837494736739572135525742023001037992695260175335123517729291);
				   Scalar b30_25(-2.11046308125864932128717300046622750300375054278936987850718);
				   Scalar b30_26(0);
				   Scalar b30_27(0.584129605671551340432988730158480872095335329645227595707052);
				   Scalar b30_28(0.711915811585189227887648262043794387578291882406745570495765);
				   Scalar b30_29(0.550499525310802324138388507020508177411414311000037561712836);
				   system->GetCurrDerivatives(k30, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b30_0 * k1[coordIndex]
						   + b30_1 * k2[coordIndex]
						   + b30_2 * k3[coordIndex]
						   + b30_3 * k4[coordIndex]
						   + b30_4 * k5[coordIndex]
						   + b30_5 * k6[coordIndex]
						   + b30_6 * k7[coordIndex]
						   + b30_7 * k8[coordIndex]
						   + b30_8 * k9[coordIndex]
						   + b30_9 * k10[coordIndex]
						   + b30_10 * k11[coordIndex]
						   + b30_11 * k12[coordIndex]
						   + b30_12 * k13[coordIndex]
						   + b30_13 * k14[coordIndex]
						   + b30_14 * k15[coordIndex]
						   + b30_15 * k16[coordIndex]
						   + b30_16 * k17[coordIndex]
						   + b30_17 * k18[coordIndex]
						   + b30_18 * k19[coordIndex]
						   + b30_19 * k20[coordIndex]
						   + b30_20 * k21[coordIndex]
						   + b30_21 * k22[coordIndex]
						   + b30_22 * k23[coordIndex]
						   + b30_23 * k24[coordIndex]
						   + b30_24 * k25[coordIndex]
						   + b30_25 * k26[coordIndex]
						   + b30_26 * k27[coordIndex]
						   + b30_27 * k28[coordIndex]
						   + b30_28 * k29[coordIndex]
						   + b30_29 * k30[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a29, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 30:
		{
				   Scalar a30(0.669986979272772921764683785505998513938845229638460353285142);
				   Scalar b31_0(0.193333333333333333333333333333333333333333333333333333333333);
				   Scalar b31_1(0);
				   Scalar b31_2(0.220000000000000000000000000000000000000000000000000000000000);
				   Scalar b31_3(-0.0800000000000000000000000000000000000000000000000000000000000);
				   Scalar b31_4(0);
				   Scalar b31_5(0);
				   Scalar b31_6(0.109993425580724703919462404865068340845119058295846426463652);
				   Scalar b31_7(-0.254297048076270161384068506997153122141835626976703920846242);
				   Scalar b31_8(0);
				   Scalar b31_9(0.865570777116694254343770343821098281832847401233011859346737);
				   Scalar b31_10(3.32416449114093083106799552786572018336860092936986407160200);
				   Scalar b31_11(0);
				   Scalar b31_12(0);
				   Scalar b31_13(-12.0102223315977933882352385148661841260301942633996815127277);
				   Scalar b31_14(0.476601466242493239430442776862061899602963782003580209476163);
				   Scalar b31_15(-29.0243011221036390525802623213654099596251221332470910692353);
				   Scalar b31_16(0);
				   Scalar b31_17(0);
				   Scalar b31_18(0);
				   Scalar b31_19(0);
				   Scalar b31_20(0);
				   Scalar b31_21(0);
				   Scalar b31_22(0);
				   Scalar b31_23(29.0243011221036390525802623213654099596251221332470910692353);
				   Scalar b31_24(-0.476601466242493239430442776862061899602963782003580209476163);
				   Scalar b31_25(12.0102223315977933882352385148661841260301942633996815127277);
				   Scalar b31_26(0);
				   Scalar b31_27(-3.32416449114093083106799552786572018336860092936986407160200);
				   Scalar b31_28(-0.865570777116694254343770343821098281832847401233011859346737);
				   Scalar b31_29(0.254297048076270161384068506997153122141835626976703920846242);
				   Scalar b31_30(-0.109993425580724703919462404865068340845119058295846426463652);
				   system->GetCurrDerivatives(k31, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b31_0 * k1[coordIndex]
						   + b31_1 * k2[coordIndex]
						   + b31_2 * k3[coordIndex]
						   + b31_3 * k4[coordIndex]
						   + b31_4 * k5[coordIndex]
						   + b31_5 * k6[coordIndex]
						   + b31_6 * k7[coordIndex]
						   + b31_7 * k8[coordIndex]
						   + b31_8 * k9[coordIndex]
						   + b31_9 * k10[coordIndex]
						   + b31_10 * k11[coordIndex]
						   + b31_11 * k12[coordIndex]
						   + b31_12 * k13[coordIndex]
						   + b31_13 * k14[coordIndex]
						   + b31_14 * k15[coordIndex]
						   + b31_15 * k16[coordIndex]
						   + b31_16 * k17[coordIndex]
						   + b31_17 * k18[coordIndex]
						   + b31_18 * k19[coordIndex]
						   + b31_19 * k20[coordIndex]
						   + b31_20 * k21[coordIndex]
						   + b31_21 * k22[coordIndex]
						   + b31_22 * k23[coordIndex]
						   + b31_23 * k24[coordIndex]
						   + b31_24 * k25[coordIndex]
						   + b31_25 * k26[coordIndex]
						   + b31_26 * k27[coordIndex]
						   + b31_27 * k28[coordIndex]
						   + b31_28 * k29[coordIndex]
						   + b31_29 * k30[coordIndex]
						   + b31_30 * k31[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a30, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 31:
		{
				   Scalar a31(0.333333333333333333333333333333333333333333333333333333333333);
				   Scalar b32_0(-0.833333333333333333333333333333333333333333333333333333333333);
				   Scalar b32_1(1.38888888888888888888888888888888888888888888888888888888889);
				   Scalar b32_2(0);
				   Scalar b32_3(0);
				   Scalar b32_4(-0.750000000000000000000000000000000000000000000000000000000000);
				   Scalar b32_5(0);
				   Scalar b32_6(-0.492529543718026304422682049114021320200214681580657784719074);
				   Scalar b32_7(0);
				   Scalar b32_8(0);
				   Scalar b32_9(0);
				   Scalar b32_10(0);
				   Scalar b32_11(0);
				   Scalar b32_12(0);
				   Scalar b32_13(0);
				   Scalar b32_14(0);
				   Scalar b32_15(0);
				   Scalar b32_16(0);
				   Scalar b32_17(0);
				   Scalar b32_18(0);
				   Scalar b32_19(0);
				   Scalar b32_20(0);
				   Scalar b32_21(0);
				   Scalar b32_22(0);
				   Scalar b32_23(0);
				   Scalar b32_24(0);
				   Scalar b32_25(0);
				   Scalar b32_26(0);
				   Scalar b32_27(0);
				   Scalar b32_28(0);
				   Scalar b32_29(0);
				   Scalar b32_30(0.492529543718026304422682049114021320200214681580657784719074);
				   Scalar b32_31(0.750000000000000000000000000000000000000000000000000000000000);
				   system->GetCurrDerivatives(k32, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b32_0 * k1[coordIndex]
						   + b32_1 * k2[coordIndex]
						   + b32_2 * k3[coordIndex]
						   + b32_3 * k4[coordIndex]
						   + b32_4 * k5[coordIndex]
						   + b32_5 * k6[coordIndex]
						   + b32_6 * k7[coordIndex]
						   + b32_7 * k8[coordIndex]
						   + b32_8 * k9[coordIndex]
						   + b32_9 * k10[coordIndex]
						   + b32_10 * k11[coordIndex]
						   + b32_11 * k12[coordIndex]
						   + b32_12 * k13[coordIndex]
						   + b32_13 * k14[coordIndex]
						   + b32_14 * k15[coordIndex]
						   + b32_15 * k16[coordIndex]
						   + b32_16 * k17[coordIndex]
						   + b32_17 * k18[coordIndex]
						   + b32_18 * k19[coordIndex]
						   + b32_19 * k20[coordIndex]
						   + b32_20 * k21[coordIndex]
						   + b32_21 * k22[coordIndex]
						   + b32_22 * k23[coordIndex]
						   + b32_23 * k24[coordIndex]
						   + b32_24 * k25[coordIndex]
						   + b32_25 * k26[coordIndex]
						   + b32_26 * k27[coordIndex]
						   + b32_27 * k28[coordIndex]
						   + b32_28 * k29[coordIndex]
						   + b32_29 * k30[coordIndex]
						   + b32_30 * k31[coordIndex]
						   + b32_31 * k32[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a31, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 32:
		{
				   Scalar a32(0.555555555555555555555555555555555555555555555555555555555556);
				   Scalar b33_0(0.111111111111111111111111111111111111111111111111111111111111);
				   Scalar b33_1(0);
				   Scalar b33_2(-0.222222222222222222222222222222222222222222222222222222222222);
				   Scalar b33_3(0);
				   Scalar b33_4(0);
				   Scalar b33_5(0);
				   Scalar b33_6(0);
				   Scalar b33_7(0);
				   Scalar b33_8(0);
				   Scalar b33_9(0);
				   Scalar b33_10(0);
				   Scalar b33_11(0);
				   Scalar b33_12(0);
				   Scalar b33_13(0);
				   Scalar b33_14(0);
				   Scalar b33_15(0);
				   Scalar b33_16(0);
				   Scalar b33_17(0);
				   Scalar b33_18(0);
				   Scalar b33_19(0);
				   Scalar b33_20(0);
				   Scalar b33_21(0);
				   Scalar b33_22(0);
				   Scalar b33_23(0);
				   Scalar b33_24(0);
				   Scalar b33_25(0);
				   Scalar b33_26(0);
				   Scalar b33_27(0);
				   Scalar b33_28(0);
				   Scalar b33_29(0);
				   Scalar b33_30(0);
				   Scalar b33_31(0);
				   Scalar b33_32(0.222222222222222222222222222222222222222222222222222222222222);
				   system->GetCurrDerivatives(k33, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b33_0 * k1[coordIndex]
						   + b33_1 * k2[coordIndex]
						   + b33_2 * k3[coordIndex]
						   + b33_3 * k4[coordIndex]
						   + b33_4 * k5[coordIndex]
						   + b33_5 * k6[coordIndex]
						   + b33_6 * k7[coordIndex]
						   + b33_7 * k8[coordIndex]
						   + b33_8 * k9[coordIndex]
						   + b33_9 * k10[coordIndex]
						   + b33_10 * k11[coordIndex]
						   + b33_11 * k12[coordIndex]
						   + b33_12 * k13[coordIndex]
						   + b33_13 * k14[coordIndex]
						   + b33_14 * k15[coordIndex]
						   + b33_15 * k16[coordIndex]
						   + b33_16 * k17[coordIndex]
						   + b33_17 * k18[coordIndex]
						   + b33_18 * k19[coordIndex]
						   + b33_19 * k20[coordIndex]
						   + b33_20 * k21[coordIndex]
						   + b33_21 * k22[coordIndex]
						   + b33_22 * k23[coordIndex]
						   + b33_23 * k24[coordIndex]
						   + b33_24 * k25[coordIndex]
						   + b33_25 * k26[coordIndex]
						   + b33_26 * k27[coordIndex]
						   + b33_27 * k28[coordIndex]
						   + b33_28 * k29[coordIndex]
						   + b33_29 * k30[coordIndex]
						   + b33_30 * k31[coordIndex]
						   + b33_31 * k32[coordIndex]
						   + b33_32 * k33[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a32, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 33:
		{
				   Scalar a33(0.111111111111111111111111111111111111111111111111111111111111);
				   Scalar b34_0(0.285835140388971558796088842163836414852927537894596466840753);
				   Scalar b34_1(0.291666666666666666666666666666666666666666666666666666666667);
				   Scalar b34_2(0.218750000000000000000000000000000000000000000000000000000000);
				   Scalar b34_3(0);
				   Scalar b34_4(0.164062500000000000000000000000000000000000000000000000000000);
				   Scalar b34_5(0);
				   Scalar b34_6(0.218194354945556658327188241581352107093288824322187941141516);
				   Scalar b34_7(0.180392898478697766863635221946775437719620053641849228562435);
				   Scalar b34_8(0);
				   Scalar b34_9(0.205713839404845018859120755122929542277570094982808905393991);
				   Scalar b34_10(0.242715791581770239970282927959446515762745971386670541948576);
				   Scalar b34_11(0.246465780813629305833609291181891407799228103869305705137021);
				   Scalar b34_12(-3.44991940790890824979834154601622662060370460614931644223924);
				   Scalar b34_13(0.228875562160036081760729060738458584294220372552740218459295);
				   Scalar b34_14(0.283290599702151415321527419056733335978436595493855789831434);
				   Scalar b34_15(3.21085125837766640960131490544236787005557320332238705967955);
				   Scalar b34_16(-0.223538777364845699920233756214162507964125230083674032084065);
				   Scalar b34_17(-0.707121157204419073518727286207487212130091231955206160635271);
				   Scalar b34_18(3.21123345150287080408174729202856500893260034443022374267639);
				   Scalar b34_19(1.40954348309669766030414474301123175769045945573548986335553);
				   Scalar b34_20(-0.151362053443742613121602276742518111090963026203676055891793);
				   Scalar b34_21(0.372350574527014276454724080214619984397121028202148298716575);
				   Scalar b34_22(0.252978746406361336722199907762141285915775728129414319261111);
				   Scalar b34_23(-3.21085125837766640960131490544236787005557320332238705967955);
				   Scalar b34_24(-0.283290599702151415321527419056733335978436595493855789831434);
				   Scalar b34_25(-0.228875562160036081760729060738458584294220372552740218459295);
				   Scalar b34_26(-0.246465780813629305833609291181891407799228103869305705137021);
				   Scalar b34_27(-0.242715791581770239970282927959446515762745971386670541948576);
				   Scalar b34_28(-0.205713839404845018859120755122929542277570094982808905393991);
				   Scalar b34_29(-0.180392898478697766863635221946775437719620053641849228562435);
				   Scalar b34_30(-0.218194354945556658327188241581352107093288824322187941141516);
				   Scalar b34_31(-0.164062500000000000000000000000000000000000000000000000000000);
				   Scalar b34_32(-0.218750000000000000000000000000000000000000000000000000000000);
				   Scalar b34_33(-0.291666666666666666666666666666666666666666666666666666666667);
				   system->GetCurrDerivatives(k34, this->globalStepIndex, this->hierarchyPhase);
#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b34_0 * k1[coordIndex]
						   + b34_1 * k2[coordIndex]
						   + b34_2 * k3[coordIndex]
						   + b34_3 * k4[coordIndex]
						   + b34_4 * k5[coordIndex]
						   + b34_5 * k6[coordIndex]
						   + b34_6 * k7[coordIndex]
						   + b34_7 * k8[coordIndex]
						   + b34_8 * k9[coordIndex]
						   + b34_9 * k10[coordIndex]
						   + b34_10 * k11[coordIndex]
						   + b34_11 * k12[coordIndex]
						   + b34_12 * k13[coordIndex]
						   + b34_13 * k14[coordIndex]
						   + b34_14 * k15[coordIndex]
						   + b34_15 * k16[coordIndex]
						   + b34_16 * k17[coordIndex]
						   + b34_17 * k18[coordIndex]
						   + b34_18 * k19[coordIndex]
						   + b34_19 * k20[coordIndex]
						   + b34_20 * k21[coordIndex]
						   + b34_21 * k22[coordIndex]
						   + b34_22 * k23[coordIndex]
						   + b34_23 * k24[coordIndex]
						   + b34_24 * k25[coordIndex]
						   + b34_25 * k26[coordIndex]
						   + b34_26 * k27[coordIndex]
						   + b34_27 * k28[coordIndex]
						   + b34_28 * k29[coordIndex]
						   + b34_29 * k30[coordIndex]
						   + b34_30 * k31[coordIndex]
						   + b34_31 * k32[coordIndex]
						   + b34_32 * k33[coordIndex]
						   + b34_33 * k34[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep *a33, probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;
		case 34:
{

	Scalar a33(0.111111111111111111111111111111111111111111111111111111111111);
	Scalar a34(1.00000000000000000000000000000000000000000000000000000000000);
	Scalar c0(0.0178571428571428571428571428571428571428571428571428571428571);
	Scalar c1(0.00585937500000000000000000000000000000000000000000000000000000);
	Scalar c2(0.0117187500000000000000000000000000000000000000000000000000000);
	Scalar c3(0);
	Scalar c4(0.0175781250000000000000000000000000000000000000000000000000000);
	Scalar c5(0);
	Scalar c6(0.0234375000000000000000000000000000000000000000000000000000000);
	Scalar c7(0.0292968750000000000000000000000000000000000000000000000000000);
	Scalar c8(0);
	Scalar c9(0.0351562500000000000000000000000000000000000000000000000000000);
	Scalar c10(0.0410156250000000000000000000000000000000000000000000000000000);
	Scalar c11(0.0468750000000000000000000000000000000000000000000000000000000);
	Scalar c12(0);
	Scalar c13(0.0527343750000000000000000000000000000000000000000000000000000);
	Scalar c14(0.0585937500000000000000000000000000000000000000000000000000000);
	Scalar c15(0.0644531250000000000000000000000000000000000000000000000000000);
	Scalar c16(0);
	Scalar c17(0.105352113571753019691496032887878162227673083080523884041670);
	Scalar c18(0.170561346241752182382120338553874085887555487802790804737501);
	Scalar c19(0.206229397329351940783526485701104894741914286259542454077972);
	Scalar c20(0.206229397329351940783526485701104894741914286259542454077972);
	Scalar c21(0.170561346241752182382120338553874085887555487802790804737501);
	Scalar c22(0.105352113571753019691496032887878162227673083080523884041670);
	Scalar c23(-0.0644531250000000000000000000000000000000000000000000000000000);
	Scalar c24(-0.0585937500000000000000000000000000000000000000000000000000000);
	Scalar c25(-0.0527343750000000000000000000000000000000000000000000000000000);
	Scalar c26(-0.0468750000000000000000000000000000000000000000000000000000000);
	Scalar c27(-0.0410156250000000000000000000000000000000000000000000000000000);
	Scalar c28(-0.0351562500000000000000000000000000000000000000000000000000000);
	Scalar c29(-0.0292968750000000000000000000000000000000000000000000000000000);
	Scalar c30(-0.0234375000000000000000000000000000000000000000000000000000000);
	Scalar c31(-0.0175781250000000000000000000000000000000000000000000000000000);
	Scalar c32(-0.0117187500000000000000000000000000000000000000000000000000000);
	Scalar c33(-0.00585937500000000000000000000000000000000000000000000000000000);
	Scalar c34(0.0178571428571428571428571428571428571428571428571428571428571);

	Scalar b34_0(0.285835140388971558796088842163836414852927537894596466840753);
	Scalar b34_1(0.291666666666666666666666666666666666666666666666666666666667);
	Scalar b34_2(0.218750000000000000000000000000000000000000000000000000000000);
	Scalar b34_3(0);
	Scalar b34_4(0.164062500000000000000000000000000000000000000000000000000000);
	Scalar b34_5(0);
	Scalar b34_6(0.218194354945556658327188241581352107093288824322187941141516);
	Scalar b34_7(0.180392898478697766863635221946775437719620053641849228562435);
	Scalar b34_8(0);
	Scalar b34_9(0.205713839404845018859120755122929542277570094982808905393991);
	Scalar b34_10(0.242715791581770239970282927959446515762745971386670541948576);
	Scalar b34_11(0.246465780813629305833609291181891407799228103869305705137021);
	Scalar b34_12(-3.44991940790890824979834154601622662060370460614931644223924);
	Scalar b34_13(0.228875562160036081760729060738458584294220372552740218459295);
	Scalar b34_14(0.283290599702151415321527419056733335978436595493855789831434);
	Scalar b34_15(3.21085125837766640960131490544236787005557320332238705967955);
	Scalar b34_16(-0.223538777364845699920233756214162507964125230083674032084065);
	Scalar b34_17(-0.707121157204419073518727286207487212130091231955206160635271);
	Scalar b34_18(3.21123345150287080408174729202856500893260034443022374267639);
	Scalar b34_19(1.40954348309669766030414474301123175769045945573548986335553);
	Scalar b34_20(-0.151362053443742613121602276742518111090963026203676055891793);
	Scalar b34_21(0.372350574527014276454724080214619984397121028202148298716575);
	Scalar b34_22(0.252978746406361336722199907762141285915775728129414319261111);
	Scalar b34_23(-3.21085125837766640960131490544236787005557320332238705967955);
	Scalar b34_24(-0.283290599702151415321527419056733335978436595493855789831434);
	Scalar b34_25(-0.228875562160036081760729060738458584294220372552740218459295);
	Scalar b34_26(-0.246465780813629305833609291181891407799228103869305705137021);
	Scalar b34_27(-0.242715791581770239970282927959446515762745971386670541948576);
	Scalar b34_28(-0.205713839404845018859120755122929542277570094982808905393991);
	Scalar b34_29(-0.180392898478697766863635221946775437719620053641849228562435);
	Scalar b34_30(-0.218194354945556658327188241581352107093288824322187941141516);
	Scalar b34_31(-0.164062500000000000000000000000000000000000000000000000000000);
	Scalar b34_32(-0.218750000000000000000000000000000000000000000000000000000000);
	Scalar b34_33(-0.291666666666666666666666666666666666666666666666666666666667);

	system->GetCurrDerivatives(k35, this->globalStepIndex, this->hierarchyPhase);
	#pragma omp parallel for
	for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
	{
		nextCoords1[coordIndex] = currCoords[coordIndex] + this->timeStep * (
			b34_0 * k1[coordIndex]
			+ b34_1 * k2[coordIndex]
			+ b34_2 * k3[coordIndex]
			+ b34_3 * k4[coordIndex]
			+ b34_4 * k5[coordIndex]
			+ b34_5 * k6[coordIndex]
			+ b34_6 * k7[coordIndex]
			+ b34_7 * k8[coordIndex]
			+ b34_8 * k9[coordIndex]
			+ b34_9 * k10[coordIndex]
			+ b34_10 * k11[coordIndex]
			+ b34_11 * k12[coordIndex]
			+ b34_12 * k13[coordIndex]
			+ b34_13 * k14[coordIndex]
			+ b34_14 * k15[coordIndex]
			+ b34_15 * k16[coordIndex]
			+ b34_16 * k17[coordIndex]
			+ b34_17 * k18[coordIndex]
			+ b34_18 * k19[coordIndex]
			+ b34_19 * k20[coordIndex]
			+ b34_20 * k21[coordIndex]
			+ b34_21 * k22[coordIndex]
			+ b34_22 * k23[coordIndex]
			+ b34_23 * k24[coordIndex]
			+ b34_24 * k25[coordIndex]
			+ b34_25 * k26[coordIndex]
			+ b34_26 * k27[coordIndex]
			+ b34_27 * k28[coordIndex]
			+ b34_28 * k29[coordIndex]
			+ b34_29 * k30[coordIndex]
			+ b34_30 * k31[coordIndex]
			+ b34_31 * k32[coordIndex]
			+ b34_32 * k33[coordIndex]
			+ b34_33 * k34[coordIndex]
		);
		nextCoords2[coordIndex] = currCoords[coordIndex] + this->timeStep * (
		c0  * k1[coordIndex]
		+ c1  * k2[coordIndex]
		+ c2  * k3[coordIndex]
		+ c3  * k4[coordIndex]
		+ c4  * k5[coordIndex]
		+ c5  * k6[coordIndex]
		+ c6  * k7[coordIndex]
		+ c7  * k8[coordIndex]
		+ c8  * k9[coordIndex]
		+ c9  * k10[coordIndex]
		+ c10  * k11[coordIndex]
		+ c11  * k12[coordIndex]
		+ c12  * k13[coordIndex]
		+ c13  * k14[coordIndex]
		+ c14  * k15[coordIndex]
		+ c15  * k16[coordIndex]
		+ c16  * k17[coordIndex]
		
		+ c17 * k18[coordIndex]
		+ c18  * k19[coordIndex]
		+ c19 * k20[coordIndex]
		+ c20  * k21[coordIndex]
		+ c21  * k22[coordIndex]
		+ c22  * k23[coordIndex]
		+ c23  * k24[coordIndex]
		+ c24  * k25[coordIndex]
		+ c25  * k26[coordIndex]
		+ c26 * k27[coordIndex]
		+ c27  * k28[coordIndex]
		+ c28  * k29[coordIndex]
		+ c29  * k30[coordIndex]
		+ c30  * k31[coordIndex]
		+ c31  * k32[coordIndex]
		+ c32  * k33[coordIndex]
		+ c33  * k34[coordIndex]
		+ c34  * k35[coordIndex]
		);
	}
	 stepError = (system->GetErrorValue(this->currTime, nextCoords1, nextCoords2, this->globalStepIndex, this->hierarchyPhase) / this->tolerance) / this->timeStep;
	 predictedStep = this->timeStep * Scalar( pow(a34 / stepError, a33));
} break;
}
}

	void AdvanceStep()
	{
		if (this->hierarchyPhase == 1)
		{
			this->currTime += this->timeStep;
		}
		system->SetCurrCoords(this->currTime, nextCoords2, this->globalStepIndex, this->hierarchyPhase);
	}

	void RevertStep()
	{
		system->SetCurrCoords(this->currTime, currCoords, this->globalStepIndex, this->hierarchyPhase);
	}

	Scalar GetLastStepError()
	{
		return stepError;
	}

	Scalar GetTimeStepPrediction()
	{
		return predictedStep;
	}

private:
	Scalar* currCoords;
	Scalar* nextCoords1;
	Scalar* nextCoords2;
	Scalar* probeCoords;
	Scalar* derivatives;

	Scalar* k1;
	Scalar* k2;
	Scalar* k3;
	Scalar* k4;
	Scalar* k5;
	Scalar* k6;
	Scalar* k7;
	Scalar *k8;
	Scalar *k9;
	Scalar *k10;
	Scalar *k11;
	Scalar *k12;
	Scalar *k13;
	Scalar *k14;
	Scalar *k15;
	Scalar *k16;
	Scalar *k17;
	Scalar *k18;
	Scalar *k19;
	Scalar *k20;
	Scalar *k21;
	Scalar *k22;
	Scalar *k23;
	Scalar *k24;
	Scalar *k25;
	Scalar *k26;
	Scalar *k27;
	Scalar *k28;
	Scalar *k29;
	Scalar *k30;
	Scalar *k31;
	Scalar *k32;
	Scalar *k33;
	Scalar *k34;
	Scalar *k35;

	Scalar stepError;
	Scalar predictedStep;

	DifferentialSystem<Scalar>* system;

};
