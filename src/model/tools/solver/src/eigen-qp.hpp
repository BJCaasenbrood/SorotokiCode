/*
    EigenQP: Fast quadradic programming template library based on Eigen.
 
    From https://github.com/jarredbarber/eigen-QP
 
    MIT License

    Copyright (c) 2017 Jarred Barber

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
#ifndef _EIGEN_QP_H_
#define _EIGEN_QP_H_

#include "Eigen/Dense"

/**
 * Solves quadradic programs with equality constraints
 *  using direct matrix factorization of the KKT system.
 */


namespace EigenQP
{

// Default tolerance levels specialized on types
template<typename t> t defTol();
template<>
    inline double defTol<double>() { return 1E-9; }
template<>
    inline float defTol<float>()   { return 1E-4f; }

/*
 * Solver for equality constrained problems.
 *  The KKT conditions are linear here, so we just
 *  invert with an LDLT decomposition.
 */
template<typename Scalar, int NVars=-1, int NEq=-1>
class QPEqSolver
{
private:
    const int n;
    const int m;

    static constexpr int NWork = 
           ((NVars == -1) || (NEq == -1)) ? -1 : (NVars+NEq);
    Eigen::Matrix<Scalar,NWork,NWork> Z;
    Eigen::Matrix<Scalar,NWork,1> C;

public:
    QPEqSolver(int n_vars=NVars, int n_const=NEq) : n(n_vars),m(n_const),
            Z(n+m,n+m), C(n+m,1)
        {
            Z.block(n,n,m,m).setZero();
        }
    void solve(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c, 
              Eigen::Matrix<Scalar,NEq,NVars> &A, Eigen::Matrix<Scalar,NEq,1> &b,
              Eigen::Matrix<Scalar,NVars,1> &x)
    {
        Z.block(0,0,n,n) = Q;
        Z.block(0,n,n,m) = A.adjoint();
        Z.block(n,0,m,n) = A;

        C.head(n) = -c;
        C.tail(m) = b;

        x = Z.ldlt().solve(C).head(n);
    }
};

/*
 * Solver for inequality constrained problems
 */
template<typename Scalar, int NVars=-1, int NIneq=-1>
class QPIneqSolver
{
    typedef Eigen::Matrix<Scalar,NVars,1> PVec;
    typedef Eigen::Matrix<Scalar,NIneq,1> DVec; // Dual (i.e., Lagrange multiplier) vector
    typedef Eigen::Matrix<Scalar,NVars,NVars> PMat;
private:
    // Problem size
    const int n;
    const int m;
    
    // Work buffers
    DVec s;
    DVec z;

    PVec rd;
    DVec rp;
    DVec rs;

    PVec dx;
    DVec ds;
    DVec dz;

    PVec x;
    
public:
    // Parameters
    Scalar tolerance;
    int    max_iters;

    QPIneqSolver(int n_vars=NVars, int n_const=NIneq) : n(n_vars),m(n_const), s(m), z(m), rd(n), rp(m), rs(m), dx(n), ds(m), dz(m)
        {
            tolerance = defTol<Scalar>();
            max_iters = 250;
            if (NVars == -1) {
                x.resize(n_vars);
            }
            if (NIneq == -1) {
                s.resize(n_const);
                z.resize(n_const);
            }
        }

    ~QPIneqSolver() {}

    void solve(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c, 
              Eigen::Matrix<Scalar,NIneq,NVars> &A, Eigen::Matrix<Scalar,NIneq,1> &b,
              Eigen::Matrix<Scalar,NVars,1> &x_out)
    {
        const Scalar eta(0.95);
        const Scalar eps = tolerance;

        // Initialization
        s.setOnes();
        z.setOnes();
        x.setZero();

        // Initial residuals. Uses fact that x=0 here.
        rd = c - A.adjoint()*z;
        rp = s + b;
        rs = (s.array()*z.array());

        const Scalar ms = Scalar(1.0)/(Scalar)m;
        Scalar mu = (Scalar)n*ms; // Initial mu based on knowing that s,z are ones.
        Scalar alpha;

        for (int iter=0; iter < max_iters; iter++)
        {
            // Precompute decompositions for this iteration
            Eigen::LLT<PMat> Gbar = (Q + A.adjoint()*((z.array()/s.array()).matrix().asDiagonal())*A).llt();

            for (int ii=0; ii < 2; ii++)
            {   
                // Prediction/correction step
                {
                    auto tmp = (rs.array() - z.array()*rp.array())/s.array();
                    dx = -Gbar.solve(rd + A.adjoint()*tmp.matrix());
                    ds = A*dx - rp;
                    dz.array() = -(rs.array() - z.array()*ds.array())/s.array();
                }

                // Compute alph,mu 
                alpha = 1.0;
                for (int jj=0; jj < m; jj++)
                {
                    Scalar a = -z(jj)/dz(jj);
                    alpha    = (a < alpha) && (a > 0) ? a : alpha;
                    a        = -s(jj)/ds(jj);
                    alpha    = (a < alpha) && (a > 0) ? a : alpha;
                }

                if (ii)
                    break; // Don't need to compute any more

                // Centering    
                Scalar mu_aff = (s + alpha*ds).dot(z+alpha*dz)*ms;
                Scalar sigma  = (mu_aff/mu); sigma *= sigma*sigma;

                // Corrector residual
                rs.array() += ds.array()*dz.array() - sigma*mu;
            }

            // Step
            alpha *= eta; // rescale step size
            x += alpha*dx;
            s += alpha*ds;
            z += alpha*dz;

            // Update residuals
            rd = Q*x + c - A.adjoint()*z;
            rp = s - A*x + b;
            rs = (s.array()*z.array());

            mu = s.dot(z)*ms;

            // Convergence test
            if ( (mu < eps) && 
                 (rd.norm() < eps) && 
                 (rs.norm() < eps) )
            {
                break;
            }
        }
        x_out = x;
    }
    public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#if 0
/**
 * General QPs with both equality and inequality constraints.
 *  This doesn't currently work.
 */
template<typename Scalar, int NVars=-1, int NEq=-1, int NIneq=-1>
class QPGenSolver
{
    // Static size for work matrix.
    static constexpr int NWork = 
        ((NVars == -1) || (NEq == -1) || (NIneq==-1)) ? -1 : (NVars+NEq+2*NIneq);
    typedef Eigen::Matrix<Scalar,NVars,1> PVec;
    typedef Eigen::Matrix<Scalar,NIneq,1> DVec; // Dual (i.e., Lagrange multiplier) vector
    typedef Eigen::Matrix<Scalar,NEq,1> EVec;   // Dual (i.e., Lagrange multiplier) vector for equality

    typedef Eigen::Matrix<Scalar,NWork,NWork> WorkBuf;
private:

    // Problem size
    const int n;
    const int mi;
    const int me;

    // Work buffers
    DVec s;
    DVec z;
    EVec y; 

    PVec rd;
    DVec rp;
    DVec rs;
    EVec ry;

    PVec dx;
    DVec ds;
    DVec dz;
    EVec dy;

    WorkBuf augSystem;
public:
    QPGenSolver(int n_vars=NVars, int n_const_eq=NEq, int n_const_ineq=NIneq) 
        : n(n_vars),mi(n_const_ineq),me(n_const_eq), 
          s(mi), z(mi), y(me), rd(n), rp(mi), rs(mi), 
          ry(me), dx(n), ds(mi), dz(mi), dy(me),
          augSystem(2*mi+me+n,2*mi+me+n)
        {
        }

    ~QPGenSolver() {}

    void solve(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c, 
              Eigen::Matrix<Scalar,NIneq,NVars> &A, Eigen::Matrix<Scalar,NIneq,1> &b,
              Eigen::Matrix<Scalar,NEq,NVars> &E, Eigen::Matrix<Scalar,NEq,1> &f,
              Eigen::Matrix<Scalar,NVars,1> &x)
    {
       
    }
};
#endif

template<typename Scalar, int NVars, int NIneq>
void quadprog(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c, 
              Eigen::Matrix<Scalar,NIneq,NVars> &A, Eigen::Matrix<Scalar,NIneq,1> &b,
              Eigen::Matrix<Scalar,NVars,1> &x)
{
    QPIneqSolver<Scalar,NVars,NIneq> qp(c.size(),b.size());
    qp.solve(Q,c,A,b,x);
}

} // End namespace
#endif
