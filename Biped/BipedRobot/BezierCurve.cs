using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

namespace BipedRobot
{
    //represent bezier curves in matrix form for arbitrary many controlpoints
    //written on the form y(t) = Tm*Bm*P;
    public class BezierCurve
    {
        private Matrix<double> _bezierMatrix;
        private int _numControlPoints;

        public BezierCurve(int n)
        {
            _numControlPoints = n;
            setBezierMatrix(n);
        }
        public void setBezierMatrix(int n)
        {
            Matrix<double> matrix = Matrix<double>.Build.Dense(n+1, n+1);
            for (int i = 0; i<= n; i++)
            {
                for(int j=0; j <=n; j++)
                {
                    if (((i + j)<n) && ((i + j) >= 0))
                    {
                        matrix[i, j] = binomialCoefficient(n, j) * binomialCoefficient(n - j, n - i - j) * Math.Pow(-1, n - i - j);
                    }
                    else
                    {
                        matrix[i, j] = 0;
                    }
                }
            }
            _bezierMatrix = matrix;
        }
        public double binomialCoefficient(double n, double k)
        {
            double r = 1;
            double d;
            if (k > n) return 0;
            for(d = 1;d<= k; d++)
            {
                r *= n--;
                r /= d;
            }
            return r;
        }

        public string functionToString()
        {
            string str = "";
            for (int i = 0; i<= _numControlPoints; i++)
            {
                double integer = 0;
                for (int j = 0;j<= _numControlPoints; j++)
                {
                    integer += _bezierMatrix[i, j];
                }
                str += integer.ToString() + "*" + "Pow(theta," + (_numControlPoints - i).ToString() + ")" + " ";
            }

            return str;
        }
    }
}
