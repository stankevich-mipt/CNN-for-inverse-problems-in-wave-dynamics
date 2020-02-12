#include "../DifferentialSolver.h"
#include <assert.h>

/*
Copy-past of DormandPrince method & light differences, based on http://sce.uhcl.edu/feagin/courses/rk10.pdf
results to Runge-Kutta 10th order method
with 17th stages
see also:
http://sce.uhcl.edu/rungekutta/
http://sce.uhcl.edu/rungekutta/rk108.txt
*/
template<typename Scalar>
class RungeSolver8_10 : public DifferentialSolver<Scalar>
{
public:
	RungeSolver8_10()
	{
	}
	void SetSystem(DifferentialSystem<Scalar> *system) override
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

		this->currTime = 0;
	}
	virtual ~RungeSolver8_10()
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

	}

	Scalar pow(Scalar a, Scalar b)
	{
		return exp(log(a) * b);
	}

	int GetPhasesCount() const override
	{
		return 17;
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
				  Scalar b1_0(0.100000000000000000000000000000000000000000000000000000000000);
				  system->GetCurrDerivatives(k1, this->globalStepIndex, this->hierarchyPhase);
				  #pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b1_0 * k1[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.000000000000000000000000000000000000000000000000000000000000), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		}break;

		case 1:
		{
				  Scalar b2_0(-0.915176561375291440520015019275342154318951387664369720564660);
				  Scalar b2_1(1.45453440217827322805250021715664459117622483736537873607016);
				  system->GetCurrDerivatives(k2, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b2_0 * k1[coordIndex] * this->timeStep
					  +b2_1 * k2[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.100000000000000000000000000000000000000000000000000000000000), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 2:
		{
				  Scalar b3_0(0.202259190301118170324681949205488413821477543637878380814562);
				  Scalar b3_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b3_2(0.606777570903354510974045847616465241464432630913635142443687);
				  system->GetCurrDerivatives(k3, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b3_0  * k1[coordIndex] * this->timeStep
						  + b3_1  * k2[coordIndex] * this->timeStep
						  + b3_2  * k3[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.539357840802981787532485197881302436857273449701009015505500), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		}break;

		case 3:
		{
				  Scalar b4_0(0.184024714708643575149100693471120664216774047979591417844635);
				  Scalar b4_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b4_2(0.197966831227192369068141770510388793370637287463360401555746);
				  Scalar b4_3(-0.0729547847313632629185146671595558023015011608914382961421311);
				  system->GetCurrDerivatives(k4, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b4_0 * k1[coordIndex] * this->timeStep
						  + b4_1 * k2[coordIndex] * this->timeStep
						  + b4_2 * k3[coordIndex] * this->timeStep
						  + b4_3 * k4[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.809036761204472681298727796821953655285910174551513523258250), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 4:
		{
				  Scalar b5_0(0.0879007340206681337319777094132125475918886824944548534041378);
				  Scalar b5_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b5_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b5_3(0.410459702520260645318174895920453426088035325902848695210406);
				  Scalar b5_4(0.482713753678866489204726942976896106809132737721421333413261);
				  system->GetCurrDerivatives(k5, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b5_0 * k1[coordIndex] * this->timeStep
						  + b5_1 * k2[coordIndex] * this->timeStep
						  + b5_2 * k3[coordIndex] * this->timeStep
						  + b5_3 * k4[coordIndex] * this->timeStep
						  + b5_4 * k5[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.309036761204472681298727796821953655285910174551513523258250), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		}break;

		case 5:
		{
				  Scalar b6_0(0.0859700504902460302188480225945808401411132615636600222593880);
				  Scalar b6_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b6_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b6_3(0.330885963040722183948884057658753173648240154838402033448632);
				  Scalar b6_4(0.489662957309450192844507011135898201178015478433790097210790);
				  Scalar b6_5(-0.0731856375070850736789057580558988816340355615025188195854775);
				  system->GetCurrDerivatives(k6, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + b6_0  * k1[coordIndex] * this->timeStep
						  + b6_1  * k2[coordIndex] * this->timeStep
						  + b6_2  * k3[coordIndex] * this->timeStep
						  + b6_3  * k4[coordIndex] * this->timeStep
						  + b6_4  * k5[coordIndex] * this->timeStep
						  + b6_5  * k6[coordIndex] * this->timeStep;
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.981074190219795268254879548310562080489056746118724882027805), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 6:
		{
				  Scalar b7_0(0.120930449125333720660378854927668953958938996999703678812621);
				  Scalar b7_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b7_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b7_3(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b7_4(0.260124675758295622809007617838335174368108756484693361887839);
				  Scalar b7_5(0.0325402621549091330158899334391231259332716675992700000776101);
				  Scalar b7_6(-0.0595780211817361001560122202563305121444953672762930724538856);
				  system->GetCurrDerivatives(k7, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b7_0  * k1[coordIndex]
						  + b7_1  * k2[coordIndex]
						  + b7_2  * k3[coordIndex]
						  + b7_3  * k4[coordIndex]
						  + b7_4  * k5[coordIndex]
						  + b7_5  * k6[coordIndex]
						  + b7_6  * k7[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.833333333333333333333333333333333333333333333333333333333333), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 7:
		{
				  Scalar b8_0(0.110854379580391483508936171010218441909425780168656559807038);
				  Scalar b8_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b8_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b8_3(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b8_4(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b8_5(-0.0605761488255005587620924953655516875526344415354339234619466);
				  Scalar b8_6(0.321763705601778390100898799049878904081404368603077129251110);
				  Scalar b8_7(0.510485725608063031577759012285123416744672137031752354067590);
				  system->GetCurrDerivatives(k8, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b8_0  * k1[coordIndex]
						  + b8_1  * k2[coordIndex]
						  + b8_2  * k3[coordIndex]
						  + b8_3  * k4[coordIndex]
						  + b8_4  * k5[coordIndex]
						  + b8_5  * k6[coordIndex]
						  + b8_6  * k7[coordIndex]
						  + b8_7  * k8[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.354017365856802376329264185948796742115824053807373968324184), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 8:
		{
				  Scalar b9_0(0.112054414752879004829715002761802363003717611158172229329393);
				  Scalar b9_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b9_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b9_3(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b9_4(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b9_5(-0.144942775902865915672349828340980777181668499748506838876185);
				  Scalar b9_6(-0.333269719096256706589705211415746871709467423992115497968724);
				  Scalar b9_7(0.499269229556880061353316843969978567860276816592673201240332);
				  Scalar b9_8(0.509504608929686104236098690045386253986643232352989602185060);
				  system->GetCurrDerivatives(k9, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b9_0  * k1[coordIndex]
						  + b9_1  * k2[coordIndex]
						  + b9_2  * k3[coordIndex]
						  + b9_3  * k4[coordIndex]
						  + b9_4  * k5[coordIndex]
						  + b9_5  * k6[coordIndex]
						  + b9_6  * k7[coordIndex]
						  + b9_7  * k8[coordIndex]
						  + b9_8  * k9[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.882527661964732346425501486979669075182867844268052119663791), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 9:
		{
				  Scalar b10_0(0.113976783964185986138004186736901163890724752541486831640341);
				  Scalar b10_1(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b10_2(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b10_3(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b10_4(0.000000000000000000000000000000000000000000000000000000000000);
				  Scalar b10_5(-0.0768813364203356938586214289120895270821349023390922987406384);
				  Scalar b10_6(0.239527360324390649107711455271882373019741311201004119339563);
				  Scalar b10_7(0.397774662368094639047830462488952104564716416343454639902613);
				  Scalar b10_8(0.0107558956873607455550609147441477450257136782823280838547024);
				  Scalar b10_9(-0.327769124164018874147061087350233395378262992392394071906457);
				  system->GetCurrDerivatives(k10, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				  for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				  {
					  probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						  b10_0  * k1[coordIndex]
						  + b10_1  * k2[coordIndex]
						  + b10_2  * k3[coordIndex]
						  + b10_3  * k4[coordIndex]
						  + b10_4  * k5[coordIndex]
						  + b10_5  * k6[coordIndex]
						  + b10_6  * k7[coordIndex]
						  + b10_7  * k8[coordIndex]
						  + b10_8  * k9[coordIndex]
						  + b10_9  * k10[coordIndex]
						  );
				  }
				  system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.642615758240322548157075497020439535959501736363212695909875), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 10:
		{
				   Scalar b11_0(0.0798314528280196046351426864486400322758737630423413945356284);
				   Scalar b11_1(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b11_2(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b11_3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b11_4(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b11_5(-0.0520329686800603076514949887612959068721311443881683526937298);
				   Scalar b11_6(-0.0576954146168548881732784355283433509066159287152968723021864);
				   Scalar b11_7(0.194781915712104164976306262147382871156142921354409364738090);
				   Scalar b11_8(0.145384923188325069727524825977071194859203467568236523866582);
				   Scalar b11_9(-0.0782942710351670777553986729725692447252077047239160551335016);
				   Scalar b11_10(-0.114503299361098912184303164290554670970133218405658122674674);
				   system->GetCurrDerivatives(k11, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b11_0  * k1[coordIndex]
						   + b11_1  * k2[coordIndex]
						   + b11_2  * k3[coordIndex]
						   + b11_3  * k4[coordIndex]
						   + b11_4  * k5[coordIndex]
						   + b11_5  * k6[coordIndex]
						   + b11_6  * k7[coordIndex]
						   + b11_7  * k8[coordIndex]
						   + b11_8  * k9[coordIndex]
						   + b11_9  * k10[coordIndex]
						   + b11_10  * k11[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.357384241759677451842924502979560464040498263636787304090125), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 11:
		{
				   Scalar b12_0(0.985115610164857280120041500306517278413646677314195559520529);
				   Scalar b12_1(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b12_2(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b12_3(0.330885963040722183948884057658753173648240154838402033448632);
				   Scalar b12_4(0.489662957309450192844507011135898201178015478433790097210790);
				   Scalar b12_5(-1.37896486574843567582112720930751902353904327148559471526397);
				   Scalar b12_6(-0.861164195027635666673916999665534573351026060987427093314412);
				   Scalar b12_7(5.78428813637537220022999785486578436006872789689499172601856);
				   Scalar b12_8(3.28807761985103566890460615937314805477268252903342356581925);
				   Scalar b12_9(-2.38633905093136384013422325215527866148401465975954104585807);
				   Scalar b12_10(-3.25479342483643918654589367587788726747711504674780680269911);
				   Scalar b12_11(-2.16343541686422982353954211300054820889678036420109999154887);
				   system->GetCurrDerivatives(k12, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b12_0  * k1[coordIndex]
						   + b12_1  * k2[coordIndex]
						   + b12_2  * k3[coordIndex]
						   + b12_3  * k4[coordIndex]
						   + b12_4  * k5[coordIndex]
						   + b12_5  * k6[coordIndex]
						   + b12_6  * k7[coordIndex]
						   + b12_7  * k8[coordIndex]
						   + b12_8  * k9[coordIndex]
						   + b12_9  * k10[coordIndex]
						   + b12_10  * k11[coordIndex]
						   + b12_11  * k12[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.117472338035267653574498513020330924817132155731947880336209), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 12:
		{
				   Scalar b13_0(0.895080295771632891049613132336585138148156279241561345991710);
				   Scalar b13_1(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b13_2(0.197966831227192369068141770510388793370637287463360401555746);
				   Scalar b13_3(-0.0729547847313632629185146671595558023015011608914382961421311);
				   Scalar b13_4(0.0000000000000000000000000000000000000000000000000000000000000);
				   Scalar b13_5(-0.851236239662007619739049371445966793289359722875702227166105);
				   Scalar b13_6(0.398320112318533301719718614174373643336480918103773904231856);
				   Scalar b13_7(3.63937263181035606029412920047090044132027387893977804176229);
				   Scalar b13_8(1.54822877039830322365301663075174564919981736348973496313065);
				   Scalar b13_9(-2.12221714704053716026062427460427261025318461146260124401561);
				   Scalar b13_10(-1.58350398545326172713384349625753212757269188934434237975291);
				   Scalar b13_11(-1.71561608285936264922031819751349098912615880827551992973034);
				   Scalar b13_12(-0.0244036405750127452135415444412216875465593598370910566069132);
				   system->GetCurrDerivatives(k13, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b13_0  * k1[coordIndex]
						   + b13_1  * k2[coordIndex]
						   + b13_2  * k3[coordIndex]
						   + b13_3  * k4[coordIndex]
						   + b13_4  * k5[coordIndex]
						   + b13_5  * k6[coordIndex]
						   + b13_6  * k7[coordIndex]
						   + b13_7  * k8[coordIndex]
						   + b13_8  * k9[coordIndex]
						   + b13_9  * k10[coordIndex]
						   + b13_10  * k11[coordIndex]
						   + b13_11  * k12[coordIndex]
						   + b13_12  * k13[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.833333333333333333333333333333333333333333333333333333333333), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 13:
		{
				   Scalar b14_0(-0.915176561375291440520015019275342154318951387664369720564660);
				   Scalar b14_1(1.45453440217827322805250021715664459117622483736537873607016);
				   Scalar b14_2(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_4(-0.777333643644968233538931228575302137803351053629547286334469);
				   Scalar b14_5(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_6(-0.0910895662155176069593203555807484200111889091770101799647985);
				   Scalar b14_7(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_8(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_9(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_10(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_11(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b14_12(0.0910895662155176069593203555807484200111889091770101799647985);
				   Scalar b14_13(0.777333643644968233538931228575302137803351053629547286334469);
				   system->GetCurrDerivatives(k14, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b14_0  * k1[coordIndex]
						   + b14_1  * k2[coordIndex]
						   + b14_2  * k3[coordIndex]
						   + b14_3  * k4[coordIndex]
						   + b14_4  * k5[coordIndex]
						   + b14_5  * k6[coordIndex]
						   + b14_6  * k7[coordIndex]
						   + b14_7  * k8[coordIndex]
						   + b14_8  * k9[coordIndex]
						   + b14_9  * k10[coordIndex]
						   + b14_10  * k11[coordIndex]
						   + b14_11  * k12[coordIndex]
						   + b14_12  * k13[coordIndex]
						   + b14_13  * k14[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.309036761204472681298727796821953655285910174551513523258250), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 14:
		{
				   Scalar b15_0(0.100000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_1(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_2(-0.157178665799771163367058998273128921867183754126709419409654);
				   Scalar b15_3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_4(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_5(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_6(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_7(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_8(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_9(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_10(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_11(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_12(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_13(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b15_14(0.157178665799771163367058998273128921867183754126709419409654);
				   system->GetCurrDerivatives(k15, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b15_0  * k1[coordIndex]
						   + b15_1  * k2[coordIndex]
						   + b15_2  * k3[coordIndex]
						   + b15_3  * k4[coordIndex]
						   + b15_4  * k5[coordIndex]
						   + b15_5  * k6[coordIndex]
						   + b15_6  * k7[coordIndex]
						   + b15_7  * k8[coordIndex]
						   + b15_8  * k9[coordIndex]
						   + b15_9  * k10[coordIndex]
						   + b15_10  * k11[coordIndex]
						   + b15_11  * k12[coordIndex]
						   + b15_12  * k13[coordIndex]
						   + b15_13  * k14[coordIndex]
						   + b15_14  * k15[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.539357840802981787532485197881302436857273449701009015505500), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 15:
		{
				   Scalar b16_0(0.181781300700095283888472062582262379650443831463199521664945);
				   Scalar b16_1(0.675000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_2(0.342758159847189839942220553413850871742338734703958919937260);
				   Scalar b16_3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_4(0.259111214548322744512977076191767379267783684543182428778156);
				   Scalar b16_5(-0.358278966717952089048961276721979397739750634673268802484271);
				   Scalar b16_6(-1.04594895940883306095050068756409905131588123172378489286080);
				   Scalar b16_7(0.930327845415626983292300564432428777137601651182965794680397);
				   Scalar b16_8(1.77950959431708102446142106794824453926275743243327790536000);
				   Scalar b16_9(0.100000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_10(-0.282547569539044081612477785222287276408489375976211189952877);
				   Scalar b16_11(-0.159327350119972549169261984373485859278031542127551931461821);
				   Scalar b16_12(-0.145515894647001510860991961081084111308650130578626404945571);
				   Scalar b16_13(-0.259111214548322744512977076191767379267783684543182428778156);
				   Scalar b16_14(-0.342758159847189839942220553413850871742338734703958919937260);
				   Scalar b16_15(-0.675000000000000000000000000000000000000000000000000000000000);
				   system->GetCurrDerivatives(k16, this->globalStepIndex, this->hierarchyPhase);
					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {
					   probeCoords[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b16_0  * k1[coordIndex]
						   + b16_1  * k2[coordIndex]
						   + b16_2  * k3[coordIndex]
						   + b16_3  * k4[coordIndex]
						   + b16_4  * k5[coordIndex]
						   + b16_5  * k6[coordIndex]
						   + b16_6  * k7[coordIndex]
						   + b16_7  * k8[coordIndex]
						   + b16_8  * k9[coordIndex]
						   + b16_9  * k10[coordIndex]
						   + b16_10  * k11[coordIndex]
						   + b16_11  * k12[coordIndex]
						   + b16_12  * k13[coordIndex]
						   + b16_13  * k14[coordIndex]
						   + b16_14  * k15[coordIndex]
						   + b16_15  * k16[coordIndex]
						   );
				   }
				   system->SetCurrCoords(this->currTime + this->timeStep * Scalar(0.100000000000000000000000000000000000000000000000000000000000), probeCoords, this->globalStepIndex, this->hierarchyPhase);
		} break;

		case 16:
		{
				   Scalar b16_0(0.181781300700095283888472062582262379650443831463199521664945);
				   Scalar b16_1(0.675000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_2(0.342758159847189839942220553413850871742338734703958919937260);
				   Scalar b16_3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_4(0.259111214548322744512977076191767379267783684543182428778156);
				   Scalar b16_5(-0.358278966717952089048961276721979397739750634673268802484271);
				   Scalar b16_6(-1.04594895940883306095050068756409905131588123172378489286080);
				   Scalar b16_7(0.930327845415626983292300564432428777137601651182965794680397);
				   Scalar b16_8(1.77950959431708102446142106794824453926275743243327790536000);
				   Scalar b16_9(0.100000000000000000000000000000000000000000000000000000000000);
				   Scalar b16_10(-0.282547569539044081612477785222287276408489375976211189952877);
				   Scalar b16_11(-0.159327350119972549169261984373485859278031542127551931461821);
				   Scalar b16_12(-0.145515894647001510860991961081084111308650130578626404945571);
				   Scalar b16_13(-0.259111214548322744512977076191767379267783684543182428778156);
				   Scalar b16_14(-0.342758159847189839942220553413850871742338734703958919937260);
				   Scalar b16_15(-0.675000000000000000000000000000000000000000000000000000000000);

				   Scalar c0(0.0333333333333333333333333333333333333333333333333333333333333);
				   Scalar c1(0.0250000000000000000000000000000000000000000000000000000000000);
				   Scalar c2(0.0333333333333333333333333333333333333333333333333333333333333);
				   Scalar c3(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar c4(0.0500000000000000000000000000000000000000000000000000000000000);
				   Scalar c5(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar c6(0.0400000000000000000000000000000000000000000000000000000000000);
				   Scalar c7(0.000000000000000000000000000000000000000000000000000000000000);
				   Scalar c8(0.189237478148923490158306404106012326238162346948625830327194);
				   Scalar c9(0.277429188517743176508360262560654340428504319718040836339472);
				   Scalar c10(0.277429188517743176508360262560654340428504319718040836339472);
				   Scalar c11(0.189237478148923490158306404106012326238162346948625830327194);
				   Scalar c12(-0.0400000000000000000000000000000000000000000000000000000000000);
				   Scalar c13(-0.0500000000000000000000000000000000000000000000000000000000000);
				   Scalar c14(-0.0333333333333333333333333333333333333333333333333333333333333);
				   Scalar c15(-0.0250000000000000000000000000000000000000000000000000000000000);
				   Scalar c16(0.0333333333333333333333333333333333333333333333333333333333333);
				   system->GetCurrDerivatives(k17, this->globalStepIndex, this->hierarchyPhase);

					#pragma omp parallel for
				   for (int coordIndex = 0; coordIndex < system->GetDimentionsCount(this->globalStepIndex, this->hierarchyPhase); coordIndex++)
				   {

					   nextCoords1[coordIndex] = currCoords[coordIndex] + this->timeStep * (
						   b16_0  * k1[coordIndex]
						   + b16_1  * k2[coordIndex]
						   + b16_2  * k3[coordIndex]
						   + b16_3  * k4[coordIndex]
						   + b16_4  * k5[coordIndex]
						   + b16_5  * k6[coordIndex]
						   + b16_6  * k7[coordIndex]
						   + b16_7  * k8[coordIndex]
						   + b16_8  * k9[coordIndex]
						   + b16_9  * k10[coordIndex]
						   + b16_10  * k11[coordIndex]
						   + b16_11  * k12[coordIndex]
						   + b16_12  * k13[coordIndex]
						   + b16_13  * k14[coordIndex]
						   + b16_14  * k15[coordIndex]
						   + b16_15  * k16[coordIndex]
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
						   );
				   }
				   stepError = (system->GetErrorValue(this->currTime, nextCoords1, nextCoords2, this->globalStepIndex, this->hierarchyPhase) / this->tolerance) / this->timeStep;
				   predictedStep = this->timeStep * Scalar(pow(Scalar(1.0) / stepError, Scalar(0.1)));
		} break;
		}
	}

	void AdvanceStep(const SolverState& solverState) override
	{
		if (this->hierarchyPhase == 1)
		{
			this->currTime += this->timeStep;
		}
		system->SetCurrCoords(this->currTime, nextCoords2, this->globalStepIndex, this->hierarchyPhase);
	}

	void RevertStep(Scalar) override 
	{
		system->SetCurrCoords(this->currTime, currCoords, this->globalStepIndex, this->hierarchyPhase);
	}

	Scalar GetLastStepError() const override
	{
		return stepError;
	}

	Scalar GetTimeStepPrediction() const override
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

	Scalar stepError;
	Scalar predictedStep;

	DifferentialSystem<Scalar>* system;
};
