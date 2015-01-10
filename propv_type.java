

public class propv_type {
	double sgc;
	int lvar;
	int mdvar;
	int klim;
	@Override
	public String toString() {
		return "propv_type [sgc=" + sgc + ", lvar=" + lvar + ", mdvar=" + mdvar
				+ ", klim=" + klim + "]";
	}
	public double getSgc() {
		return sgc;
	}
	public void setSgc(double sgc) {
		this.sgc = sgc;
	}
	public int getLvar() {
		return lvar;
	}
	public void setLvar(int lvar) {
		this.lvar = lvar;
	}
	public int getMdvar() {
		return mdvar;
	}
	public void setMdvar(int mdvar) {
		this.mdvar = mdvar;
	}
	public int getKlim() {
		return klim;
	}
	public void setKlim(int klim) {
		this.klim = klim;
	}
}
