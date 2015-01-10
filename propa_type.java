

import java.util.Arrays;

public class propa_type {
	double dlsa;
	double dx;
	double ael;
	double ak1;
	double ak2;
	double aed;
	double emd;
	double aes;
	double ems;
	double[] dls = new double[2];
	double dla;
	double tha;
	@Override
	public String toString() {
		return "propa_type [dlsa=" + dlsa + ", dx=" + dx + ", ael=" + ael
				+ ", ak1=" + ak1 + ", ak2=" + ak2 + ", aed=" + aed + ", emd="
				+ emd + ", aes=" + aes + ", ems=" + ems + ", dls="
				+ Arrays.toString(dls) + ", dla=" + dla + ", tha=" + tha + "]";
	}
	public double getDlsa() {
		return dlsa;
	}
	public void setDlsa(double dlsa) {
		this.dlsa = dlsa;
	}
	public double getDx() {
		return dx;
	}
	public void setDx(double dx) {
		this.dx = dx;
	}
	public double getAel() {
		return ael;
	}
	public void setAel(double ael) {
		this.ael = ael;
	}
	public double getAk1() {
		return ak1;
	}
	public void setAk1(double ak1) {
		this.ak1 = ak1;
	}
	public double getAk2() {
		return ak2;
	}
	public void setAk2(double ak2) {
		this.ak2 = ak2;
	}
	public double getAed() {
		return aed;
	}
	public void setAed(double aed) {
		this.aed = aed;
	}
	public double getEmd() {
		return emd;
	}
	public void setEmd(double emd) {
		this.emd = emd;
	}
	public double getAes() {
		return aes;
	}
	public void setAes(double aes) {
		this.aes = aes;
	}
	public double getEms() {
		return ems;
	}
	public void setEms(double ems) {
		this.ems = ems;
	}
	public double[] getDls() {
		return dls;
	}
	public void setDls(double[] dls) {
		this.dls = dls;
	}
	public double getDla() {
		return dla;
	}
	public void setDla(double dla) {
		this.dla = dla;
	}
	public double getTha() {
		return tha;
	}
	public void setTha(double tha) {
		this.tha = tha;
	}
}
