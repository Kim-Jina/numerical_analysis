#include <iostream>
#include <cmath>

using namespace std;

int main(void){
	short x1, x2, x3;					// three points of equation
	short X1, X2;						// variables for replacement
	short f_x1, f_x3;					// f(x1), f(x3)
	short t_value;						// tolerance value
	short gap;							// gap between two points
	short multi = 1;					// multiply
	short div1 = 1, div2 = 1;			//	divide
	short p_f_x1 = 0, p_f_x3 = 0;		// previous value of f_x1 and f_x3 to remove overflow
	
	// print equation
	cout << "y=x^3+4x^2-10=0" << endl;
	
	// input x1, x2
	cout << "input x1, x2 : ";
	cin >> x1 >> x2;

	// input tolerance value
	cout << "input tolerance value : ";
	cin >> t_value;

	gap = abs(x1 - x2);		// gab is scale between x1 and x2
	
	// find a solution
	while (0.5*gap / multi >= 1 / (float)t_value){
		multi *= 2;				// mult is multi*2 (to remove floating-point)
		X1 = multi*x1 / div1;	// replace multi*x1/div1 by X1
		X2 = multi*x2 / div2;	// replace multi*x2/div2 by X2
		x3 = (X1 + X2) / 2;		// x3 is an average value of X1 and X2
		
		// compute values of X1/multi and x3/multi in equation and multiply multi*multi*multi
		f_x1 = X1*X1*X1 + 4 * multi*X1*X1 - 10 * multi*multi*multi;
		f_x3 = x3*x3*x3 + 4 * multi*x3*x3 - 10 * multi*multi*multi;
		
		// overflow
		if (f_x1*p_f_x1 > 0){						// if previous value and present value are same sign
			if (abs(f_x1) < abs(p_f_x1))			// if the absolute value of f_x1 is less than the absolut value of p_f_x1
				break;
		}

		if (f_x3*p_f_x3>0){							// if previous value and present value are same sign
			if (abs(f_x3) < abs(p_f_x3))			// if the absolute value of f_x3 is less than the absolut value of p_f_x3
				break;
		}

		if (f_x1*f_x3 < 0){						// if f(x3) of opposite sign of f(x1)
			x2 = x3;							// x2 moves into x3
			div2 = multi;						// div2 is multi
			gap = abs(multi*x1 - x2);			// calculate and lay scale between multi*x1 and x2 in gap 
		}
		else{									// if f(x3) of same sign of f(x1)
			x1 = x3;							// x1 moves into x3
			div1 = multi;						// div1 is multi
			gap = abs(x1 - multi*x2);			// calculate and lay scale between x1 and multi*x2 in gap
		}
		
		p_f_x1 = f_x1;		// p_f_x1 is the value of f_x1
		p_f_x3 = f_x3;		// p_f_x3 is the value of f_x3
	}
	
	// print result
	cout << "result of short : " << x3 << "/" << multi << " = " << x3 / (float)multi << endl;

	return 0;
}