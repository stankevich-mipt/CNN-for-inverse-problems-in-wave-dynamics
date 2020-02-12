// Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2012

template <class Space>
class ElasticForceCorrector
{
	public:
    typedef typename Space::ElasticMaterial Node;
    typedef typename Space::Vector3 vector;
    typedef typename Space::SMatrix smatrix;
    typedef typename Space::ElasticMaterialSettings Material;
    typedef typename Space::Scalar real;

		ElasticForceCorrector()
		{
			I.identity();
		}

		void correct(Node *e, vector n, vector force) {
			const vector z = -(e->s * n) - force;
			const smatrix N00 = tsquare(n);
			e->v += -z * irhoc2 + irhoc1c2 * scalar_product(z, n) * n;
			e->s += sproduct(z, n) +
				scalar_product(z, n) / (la + 2 * mu) * (la * I - 2 * (la + mu) * N00);
		}
		
		void correct(Node *e, vector n) {
			correct(e, n, this->force);
		}

		void setForce(vector i_force) { force = i_force; }
		vector getForce() { return force; }
		void setMaterial(const Material *mat)
		{
			c1 = mat->c1();
			real c2 = mat->c2();
			rho = mat->rho();
			mu = rho * c2 * c2;
			la = rho * c1 * c1 - 2.0 * mu;
			irhoc2 = 1.0 / rho / c2;
			irhoc1c2 = 1.0 / rho * (1.0 / c2 - 1.0 / c1);
		}
	public:
		vector force;
		smatrix I;
		real irhoc2; // 1 / (rho * c2)
		real la; // la
		real mu; // mu
		real irhoc1c2; // 1 / rho * (1 / c2 - 1 / c1)
		real c1;
		real rho;
};

/*template <class Elastic, class Material, class Vector, class SMatrix>
class ElasticVelocityCorrector
{
	public:
    typedef typename Space::Elastic Node;
    typedef typename Space::Vector3 vector;
    typedef typename Space::SMatrix smatrix;

		ElasticVelocityCorrector()
		{
			I.identity();
		}

		void correct(Node *e, vector n, vector v) {
			const vector z = e->v - v;
			const smatrix N00 = tsquare(n);
			e->v = v;
			e->s += rho * (dot(z, n) * ((c1 - 2.0 * c2 - c3) * N00 + c3 * I) + c2 * sproduct(z, n));
		}
		
		void correct(Node *e, vector n) {
			correct(e, n, this->v);
		}

		void setVelocity(vector v) { this->v = v; }
		vector getVelocity() { return v; }
		void setMaterial(Material *mat)
		{
			c1 = mat->getValue("c1");
			c2 = mat->getValue("c2");
			rho = mat->getValue("rho");
			mu = rho * c2 * c2;
			la = rho * c1 * c1 - 2.0 * mu;
			irhoc2 = 1.0 / rho / c2;
			irhoc1c2 = 1.0 / rho * (1.0 / c2 - 1.0 / c1);
			c3 = c1 * la / (la + 2.0 * mu);
		}
	public:
		vector v;
		smatrix I;
		real irhoc2; // 1 / (rho * c2)
		real la; // la
		real mu; // mu
		real irhoc1c2; // 1 / rho * (1 / c2 - 1 / c1)
		real rho; // rho
		real c1, c2, c3;
};


template <class Elastic, class Vector, class SMatrix>
class ElasticMixedCorrector
{
	public:
    typedef typename Space::Elastic Node;
    typedef typename Space::Vector3 vector;
    typedef typename Space::SMatrix smatrix;

		ElasticMixedCorrector() {}

		void correct(Node *e, vector n, vector ft, real vp) {
			vector force = ft + 
			               (-fc.rho * fc.c1 * vp +
			               dot((-e->s * n - fc.rho * fc.c1 * e->v - ft), n)) * n;
			fc.correct(e, n, force);
		}
		
		void correct(Node *e, vector n) {
			correct(e, n, ft, vp);
		}

		void setForce(vector i_force) { ft = i_force; }
		vector getForce() { return ft; }
		void setVP(real i_vp) { vp = i_vp; }
		real getVP() { return vp; }
		void setMaterial(Material *mat)
		{
			fc.setMaterial(mat);
		}
	private:
		vector ft;
		real vp;
		ElasticForceCorrector<TNode> fc;
};

template <class Elastic, class Vector, class SMatrix>
class ElasticMixed2Corrector
{
	public:
    typedef typename Space::Elastic Node;
    typedef typename Space::Vector3 vector;
    typedef typename Space::SMatrix smatrix;

		ElasticMixed2Corrector() {}

		void correct(Node *e, vector n, real fp, vector vt) {
			vector v = vt +
			           n * (dot(e->v - vt, n) + 1.0 / vc.rho / vc.c1 * (dot(n, e->s * n) - fp));
			vc.correct(e, n, v);
		}
		
		void correct(Node *e, vector n) {
			correct(e, n, fp, vt);
		}

		void setFP(real i_force) { fp = i_force; }
		real getFP() { return fp; }
		void setVT(vector i_vp) { vt = i_vp; }
		vector getVT() { return vt; }
		void setMaterial(Material *mat) { vc.setMaterial(mat); }
	private:
		real fp;
		vector vt;
		ElasticVelocityCorrector<TNode> vc;
};

template <class Elastic, class Vector, class SMatrix>
class ElasticAbsorptionCorrector
{
	public:
    typedef typename Space::Elastic Node;
    typedef typename Space::Vector3 vector;
    typedef typename Space::SMatrix smatrix;

		ElasticAbsorptionCorrector()
		{
			vector n; n.null();
			fc.setForce(n);
			vc.setVelocity(n);
		}

		void correct(Node *e, vector n) {
			Node fs = *e;
			Node vs = *e;
			fc.correct(&fs, n);
			vc.correct(&vs, n);
			e->s = 0.5 * (fs.s + vs.s);
			e->v = 0.5 * (fs.v + vs.v);
		}
		void setMaterial(Material *mat)
		{
			vc.setMaterial(mat);
			fc.setMaterial(mat);
		}
	private:
		ElasticVelocityCorrector<Elastic, Vector, SMatrix> vc;
		ElasticForceCorrector<Elastic, Vector, SMatrix> fc;
};*/
