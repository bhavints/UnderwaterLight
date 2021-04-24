typedef double Doub;    // 64 - bit floating point
typedef int Int;        // 32 - bit signed integer

#define PI (float) 3.14159265358979323846

#define ASYMMETRYPARAM (float) 0.5
// FOUND HERE: https://stackoverflow.com/questions/16512817/numerical-integration-in-c

class LightFunctor
{
	private:
		double u;
		double v;
		double beta;
		double rho;

	public:
		LightFunctor(double uval, double vval, double bval, double rval) : u(uval), v(vval), beta(bval), rho(rval) {  }

		// This operator overloading enables calling
		// operator function () on objects of increment
		double operator () (Doub x) {
			Doub uprime = x;
			// Calculate phase angle function F
			double unsqDistFromLight = pow(uprime, 2) + pow(v, 2);
			double distFromLight = sqrt(unsqDistFromLight);
			double cosFromLight = -uprime / distFromLight;
			double angleFromLight = acos(cosFromLight);
			// Now calculate phase function
			// Several options, will use Henyey and Greenstein
			// https://pbr-book.org/3ed-2018/Volume_Scattering/Phase_Functions
			double phaseConstant = 1.0 / (4 * PI);
			double phaseValue = (1 - pow(ASYMMETRYPARAM, 2));
			double phaseDivisor = (1 + pow(ASYMMETRYPARAM, 2) + 2 * ASYMMETRYPARAM * cos(angleFromLight));
			phaseDivisor = pow(phaseDivisor, 1.5);
			double FinalPhase = phaseConstant * (phaseValue / phaseDivisor);
			
			double exponential = distFromLight + uprime - u;
			exponential *= beta * -1 * rho;
			exponential = exp(exponential);
			exponential /= unsqDistFromLight;

			double IntegrandVal = rho * FinalPhase * exponential;

			return IntegrandVal;
		}
};

struct Quadrature
{
	//Abstract base class for elementary quadrature algorithms.
	Int n; // Current level of refinement.

	virtual Doub next() = 0;
	//Returns the value of the integral at the nth stage of refinement. 
	//The function next() must be defined in the derived class.
};

template<class T>
struct Trapzd : Quadrature
{
	Doub a, b, s; // Limits of integration and current value of integral.
	T &func;

	Trapzd() { };

	// func is function or functor to be integrated between limits: a and b 
	Trapzd(T &funcc, const Doub aa, const Doub bb)
		: func(funcc), a(aa), b(bb)
	{
		n = 0;
	}

	// Returns the nth stage of refinement of the extended trapezoidal rule. 
	// On the first call (n = 1), the routine returns the crudest estimate  
	// of integral of f x / dx in [a,b]. Subsequent calls set n=2,3,... and
	// improve the accuracy by adding 2n - 2 additional interior points.
	Doub next()
	{
		Doub x, tnm, sum, del;
		Int it, j;
		n++;

		if (n == 1)
		{
			return (s = 0.5 * (b - a) * (func(a) + func(b)));
		}
		else
		{
			for (it = 1, j = 1; j < n - 1; j++)
			{
				it <<= 1;
			}
			tnm = it;
			// This is the spacing of the points to be added.          
			del = (b - a) / tnm;
			x = a + 0.5 * del;

			for (sum = 0.0, j = 0; j < it; j++, x += del)
			{
				sum += func(x);
			}
			// This replaces s by its refined value.  
			s = 0.5 * (s + (b - a) * sum / tnm);
			return s;
		}
	}
};

template<class T>
Doub qtrap(T &func, const Doub a, const Doub b, const Doub eps = 1.0e-10)
{
	// Returns the integral of the function or functor func from a to b. 
	// The constants EPS can be set to the desired fractional accuracy and    
	// JMAX so that 2 to the power JMAX-1 is the maximum allowed number of   
	// steps. Integration is performed by the trapezoidal rule.

	const Int JMAX = 20;
	Doub s, olds = 0.0; // Initial value of olds is arbitrary.

	Trapzd<T> t(func, a, b);

	for (Int j = 0; j < JMAX; j++)
	{
		s = t.next();

		if (j > 5) // Avoid spurious early convergence.
		{
			if (abs(s - olds) < eps * abs(olds) || (s == 0.0 && olds == 0.0))
			{
				return s;
			}
		}
		olds = s;
	}
	throw("Too many steps in routine qtrap");
}