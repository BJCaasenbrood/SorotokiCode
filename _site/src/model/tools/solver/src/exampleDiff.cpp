#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/AutoDiff>
#include <iostream>

// Computes the energy terms into costs. The first part contains termes related to edge length variation,
// and the second to angle variations.
template<typename Scalar>
void ikLikeCosts(const Eigen::Matrix<Scalar,Eigen::Dynamic,1>& curve, const Eigen::VectorXd& targetAngles, const Eigen::VectorXd& targetLengths, double beta,
                 Eigen::Matrix<Scalar,Eigen::Dynamic,1>& costs)
{
    using namespace Eigen;
    using std::atan2;

    typedef Matrix<Scalar,2,1> Vec2;
    int nb = curve.size()/2;

    costs.setZero();
    for(int k=1;k<nb-1;++k)
    {
        Vec2 pk0 = curve.template segment<2>(2*(k-1));
        Vec2 pk  = curve.template segment<2>(2*k);
        Vec2 pk1 = curve.template segment<2>(2*(k+1));

        if(k+1<nb-1) {
            costs((nb-2)+(k-1)) = (pk1-pk).norm() - targetLengths(k-1);
        }
       

        Vec2 v0 = (pk-pk0).normalized();
        Vec2 v1 = (pk1-pk).normalized();

        costs(k-1) = beta * (atan2(-v0.y() * v1.x() + v0.x() * v1.y(), v0.x() * v1.x() + v0.y() * v1.y()) - targetAngles(k-1));
    }
}

// Generic functor
template<typename _Scalar>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
                          };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    const int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; } // number of degree of freedom (= 2*nb_vertices)
    int values() const { return m_values; } // number of energy terms (= nb_vertices + nb_edges)

    // you should define that in the subclass :
    //    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

// Specialized functor warping the ikLikeCosts function
struct iklike_functor : Functor<double>
{
    typedef Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,Eigen::Dynamic,1> > ADS;
    typedef Eigen::Matrix<ADS, Eigen::Dynamic, 1> VectorXad;

    // pfirst and plast are the two extremities of the curve
    iklike_functor(const Eigen::VectorXd& targetAngles, const Eigen::VectorXd& targetLengths, double beta, const Eigen::Vector2d pfirst, const Eigen::Vector2d plast)
        :   Functor<double>(targetAngles.size()*2-4,targetAngles.size()*2-1),
            m_targetAngles(targetAngles), m_targetLengths(targetLengths), m_beta(beta),
            m_pfirst(pfirst), m_plast(plast)
    {}

    // input = x = {  ..., x_i, y_i, ....}
    // output = fvec = the value of each term
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
    {
        using namespace Eigen;
        VectorXd curves(this->inputs()+8);

        curves.segment(4,this->inputs()) = x;

        Vector2d d(1,0);
        curves.segment<2>(0)                   = m_pfirst - d;
        curves.segment<2>(2)                   = m_pfirst;
        curves.segment<2>(this->inputs()+4)    = m_plast;
        curves.segment<2>(this->inputs()+6)    = m_plast + d;

        ikLikeCosts(curves, m_targetAngles, m_targetLengths, m_beta, fvec);
        return 0;
    }

    // Compute the jacobian into fjac for the current solution x
    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac)
    {
        using namespace Eigen;
        VectorXad curves(this->inputs()+8);

        // Compute the derivatives of each degree of freedom
        // -> grad( x_i ) = (0, ..., 0, 1, 0, ..., 0) ; 1 is in position i
        for(int i=0; i<this->inputs();++i)
            curves(4+i) = ADS(x(i), this->inputs(), i);

        Vector2d d(1,0);
        curves.segment<2>(0)                   = (m_pfirst - d).cast<ADS>();
        curves.segment<2>(2)                   = (m_pfirst).cast<ADS>();
        curves.segment<2>(this->inputs()+4)    = (m_plast).cast<ADS>();
        curves.segment<2>(this->inputs()+6)    = (m_plast + d).cast<ADS>();

        VectorXad v(this->values());

        ikLikeCosts(curves, m_targetAngles, m_targetLengths, m_beta, v);

        // copy the gradient of each energy term into the Jacobian
        for(int i=0; i<this->values();++i)
            fjac.row(i) = v(i).derivatives();

        return 0;
    }

    const Eigen::VectorXd& m_targetAngles;
    const Eigen::VectorXd& m_targetLengths;
    double m_beta;
    Eigen::Vector2d m_pfirst, m_plast;
};



void draw_vecX(const Eigen::VectorXd& res)
{
  for( int i = 0; i < (res.size()/2) ; i++ )
  {
    std::cout << res.segment<2>(2*i).transpose() << "\n";
  }
}


using namespace Eigen;

int main()
{
    Eigen::Vector2d pfirst(-5., 0.);
    Eigen::Vector2d plast ( 5., 0.);

    // rest pose is a straight line starting between first and last point
    const int nb_points = 30;
    Eigen::VectorXd targetAngles (nb_points);
    targetAngles.fill(0);

    Eigen::VectorXd targetLengths(nb_points-1);
    double val = (pfirst-plast).norm() / (double)(nb_points-1);
    targetLengths.fill(val);


    // get initial solution
    Eigen::VectorXd x((nb_points-2)*2);
    for(int i = 1; i < (nb_points - 1); i++)
    {
        double s = (double)i / (double)(nb_points-1);
        x.segment<2>((i-1)*2) = plast * s + pfirst * (1. - s);
    }

    // move last point
    plast = Eigen::Vector2d(4., 1.);

    // Create the functor object
    iklike_functor func(targetAngles, targetLengths, 0.1, pfirst, plast);

    // construct the solver
    Eigen::LevenbergMarquardt<iklike_functor> lm(func);

    // adjust tolerance
    lm.parameters.ftol *= 1e-2;
    lm.parameters.xtol *= 1e-2;
    lm.parameters.maxfev = 2000;


    int a = lm.minimize(x);
    std::cerr << "info = " << a  << " " << lm.nfev << " " << lm.njev << " "  <<  "\n";

    std::cout << pfirst.transpose() << "\n";
    draw_vecX( x);
    std::cout << plast.transpose() << "\n";
}