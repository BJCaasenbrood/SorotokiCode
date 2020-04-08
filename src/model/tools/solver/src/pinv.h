


void pinv(Mxf &A){
	Mxf tmp
	JacobiSVD <MatrixXf>svd(A,ComputeThinU | ComputeThinV);
	tmp = svd.matrixV() * ( (svd.singularValues().array().abs() > 1e-2*A.norm()).select(svd.singularValues().array().inverse(), 0) ).matrix().asDiagonal() * svd.matrixU().adjoint();

}
