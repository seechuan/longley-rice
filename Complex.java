
public class Complex {
	private double re;
	private double im;

	public Complex(double real, double imag) {
		this.re = real;
		this.im = imag;
	}

	public Complex() {
	}

	public double real() {
		return re;
	}

	public double imag() {
		return im;
	}

	public double abs() {
		return Math.hypot(re, im);
	}

	public String toString() {
		if (im >= 0)
			return re + "+" + im + "i";
		else
			return re + "-" + -im + "i";
	}

	// ==============================
	// Complex number arithmetic ...
	// ==============================

	// Complex negative of Complex number ...

	public Complex Negate() {
		Complex negate = new Complex();

		negate.re = -re;
		negate.im = -im;

		return (negate);
	}

	// Compute sum of two Complex numbers cA + cB.....

	public Complex Add(Complex cB) {
		Complex sum = new Complex();

		sum.re = re + cB.re;
		sum.im = im + cB.im;

		return (sum);
	}

	// Compute difference of two Complex numbers cA - cB.....

	public Complex Sub(Complex cB) {
		Complex diff = new Complex();

		diff.re = re - cB.re;
		diff.im = im - cB.im;

		return (diff);
	}

	// Compute product of two Complex numbers cA * cB.....

	public Complex Mult(Complex cB) {
		Complex prod = new Complex();

		prod.re = re * cB.re - im * cB.im;
		prod.im = im * cB.re + re * cB.im;

		return (prod);
	}

	// Compute divisor of two Complex numbers cA / cB.....

	public Complex Div(Complex cB) {
		Complex div = new Complex();
		double dR, dDen;

		if (Math.abs(cB.re) >= Math.abs(cB.im)) {
			dR = cB.im / cB.re;
			dDen = cB.re + dR * cB.im;
			div.re = (re + dR * im) / dDen;
			div.im = (im - dR * re) / dDen;
		} else {
			dR = cB.re / cB.im;
			dDen = cB.im + dR * cB.re;
			div.re = (dR * re + im) / dDen;
			div.im = (dR * im - re) / dDen;
		}

		return (div);
	}

	// Scale Complex number by double precision number.....

	public Complex Scale(double dFactor) {
		Complex scale = new Complex();

		scale.re = dFactor * re;
		scale.im = dFactor * im;

		return (scale);
	}

	// Compute Complex number conjugate....

	public Complex Conjugate() {
		Complex conj = new Complex();

		conj.re = re;
		conj.im = -im;

		return (conj);
	}

	// Compute square root of Complex number ....

	public Complex Sqrt() {
		Complex csqrt = new Complex();
		double dX, dY, dW, dR;

		if ((re == 0) && (im == 0.0)) {
			csqrt.re = 0.0;
			csqrt.im = 0.0;
			return (csqrt);
		}

		dX = Math.abs(re);
		dY = Math.abs(im);

		if (dX >= dY) {
			dR = dY / dX;
			dW = Math.sqrt(dX)
					* Math.sqrt(0.5 * (1.0 + Math.sqrt(1 + dR * dR)));
		} else {
			dR = dX / dY;
			dW = Math.sqrt(dY) * Math.sqrt(0.5 * (dR + Math.sqrt(1 + dR * dR)));
		}

		if (re >= 0.0) {
			csqrt.re = dW;
			csqrt.im = im / (2.0 * dW);
		} else {
			csqrt.im = (im > 0.0) ? dW : -dW;
			csqrt.re = im / (2.0 * csqrt.im);
		}

		return (csqrt);
	}

	/**
	 * Argument of this Complex number (the angle in radians with the x-axis in
	 * polar coordinates).
	 * 
	 * @return arg(z) where z is this Complex number.
	 */
	public double arg() {
		return Math.atan2(im, re);
	}

	/**
	 * Complex exponential (doesn't change this Complex number).
	 * 
	 * @return exp(z) where z is this Complex number.
	 */
	public Complex exp() {
		return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re)
				* Math.sin(im));
	}

	/**
	 * Principal branch of the Complex logarithm of this Complex number.
	 * (doesn't change this Complex number). The principal branch is the branch
	 * with -pi < arg <= pi.
	 * 
	 * @return log(z) where z is this Complex number.
	 */
	public Complex log() {
		return new Complex(Math.log(this.abs()), this.arg());
	}

	// Real cosh function (used to compute Complex trig functions)
	private double cosh(double theta) {
		return (Math.exp(theta) + Math.exp(-theta)) / 2;
	}

	// Real sinh function (used to compute Complex trig functions)
	private double sinh(double theta) {
		return (Math.exp(theta) - Math.exp(-theta)) / 2;
	}

	/**
	 * Sine of this Complex number (doesn't change this Complex number). <br>
	 * sin(z) = (exp(i*z)-exp(-i*z))/(2*i).
	 * 
	 * @return sin(z) where z is this Complex number.
	 */
	public Complex sin() {
		return new Complex(cosh(im) * Math.sin(re), sinh(im) * Math.cos(re));
	}

	/**
	 * Cosine of this Complex number (doesn't change this Complex number). <br>
	 * cos(z) = (exp(i*z)+exp(-i*z))/ 2.
	 * 
	 * @return cos(z) where z is this Complex number.
	 */
	public Complex cos() {
		return new Complex(cosh(im) * Math.cos(re), -sinh(im) * Math.sin(re));
	}

	/**
	 * Hyperbolic sine of this Complex number (doesn't change this Complex
	 * number). <br>
	 * sinh(z) = (exp(z)-exp(-z))/2.
	 * 
	 * @return sinh(z) where z is this Complex number.
	 */
	public Complex sinh() {
		return new Complex(sinh(re) * Math.cos(im), cosh(re) * Math.sin(im));
	}

	/**
	 * Hyperbolic cosine of this Complex number (doesn't change this Complex
	 * number). <br>
	 * cosh(z) = (exp(z) + exp(-z)) / 2.
	 * 
	 * @return cosh(z) where z is this Complex number.
	 */
	public Complex cosh() {
		return new Complex(cosh(re) * Math.cos(im), sinh(re) * Math.sin(im));
	}

	/**
	 * Tangent of this Complex number (doesn't change this Complex number). <br>
	 * tan(z) = sin(z)/cos(z).
	 * 
	 * @return tan(z) where z is this Complex number.
	 */
	public Complex tan() {
		return (this.sin()).Div(this.cos());
	}

	public static void main(String[] args) {
		Complex c = new Complex(1, 2);
		Complex d = new Complex(2, 3);
		
		System.out.println(c.Negate().toString());
		System.out.println(c.Add(d).toString());
		System.out.println(c.Sub(d).toString());
		System.out.println(c.Mult(d).toString());
		System.out.println(c.Div(d).toString());
		System.out.println(c.Scale(Math.PI).toString());
		System.out.println(c.Conjugate().toString());
		System.out.println(c.Sqrt().toString());
		System.out.println(c.arg());
		System.out.println(c.exp());
		System.out.println(c.log());
		System.out.println(c.sin().toString());
		System.out.println(c.cos().toString());
		System.out.println(c.tan().toString());
		System.out.println(c.sinh(Math.PI));
		System.out.println(c.cosh(Math.PI));
		System.out.println(c.sinh().toString());
		System.out.println(c.cosh().toString());
		
		
		
		
	}
}
