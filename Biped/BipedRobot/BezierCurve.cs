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
using System.Globalization;

namespace BipedRobot
{
    //represent bezier curves in matrix form for arbitrary many controlpoints
    //written on the form y(t) = Tm*Bm*P;
    public class BezierCurve
    {
        private Matrix<double> _bezierMatrix;
        private int _numControlPoints;
        private Vector<double> _q1controlpoints;
        private Vector<double> _q2controlpoints;
        private Vector<double> _q3controlpoints;

        public BezierCurve(int n, BRgait gait)
        {
            _numControlPoints = n;
            setBezierMatrix(n);
            setControlPoints(gait);
        }
        public void setBezierMatrix(int n)
        {
            Matrix<double> matrix = Matrix<double>.Build.Dense(n, n);
            for (int i = 0; i< n; i++)
            {
                for(int j=0; j <n; j++)
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
        public void setControlPoints(BRgait gait)
        {
            _q1controlpoints = gait.gaitParam.gaitparameters.Item1;
            _q2controlpoints = gait.gaitParam.gaitparameters.Item2;
            _q3controlpoints = gait.gaitParam.gaitparameters.Item3;
        }
        private double binomialCoefficient(double n, double k)
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

        public Tuple<Dictionary<string,double>,string> phi1ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q1controlpoints;
            Dictionary<string, double> parameters = new Dictionary<string, double>();
            for (int i = 0; i < valueVector.Count; i++)
            {
                str += "P"+i.ToString()+" * " + "(theta^" + (_numControlPoints - 1 - i) + ")" + "+";
                parameters.Add("P" + i.ToString(), valueVector[i]);
            }
            str += "0";
            str = str.Replace("+-", "-");
            return new Tuple<Dictionary<string, double>, string>(parameters, str);
        }
        public Tuple<Dictionary<string, double>, string> phi2ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q2controlpoints;
            Dictionary<string, double> parameters = new Dictionary<string, double>();
            for (int i = 0; i < valueVector.Count; i++)
            {
                str += "P" + (i + _numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i) + ")" + "+";
                parameters.Add("P" + (i + _numControlPoints).ToString(), valueVector[i]);
            }
            str += "0";
            str = str.Replace("+-", "-");
            return new Tuple<Dictionary<string, double>, string>(parameters, str);
        }
        public Tuple<Dictionary<string, double>, string> phi3ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q3controlpoints;
            Dictionary<string, double> parameters = new Dictionary<string, double>();
            for (int i = 0; i < valueVector.Count; i++)
            {
                str += "P" + (i + 2*_numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i) + ")" + "+";
                parameters.Add("P" + (i + 2*_numControlPoints).ToString(), valueVector[i]);
            }
            str += "0";
            str = str.Replace("+-", "-");
            return new Tuple<Dictionary<string, double>, string>(parameters, str);
        }

        public string dphi1ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q1controlpoints;
            for (int i = 0; i < valueVector.Count-1; i++)
            {
                str += (_numControlPoints - 1 - i).ToString() + "*" + "P" + i.ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 1) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }
        public string dphi2ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q2controlpoints;
            for (int i = 0; i < valueVector.Count - 1; i++)
            {
                str += (_numControlPoints - 1 - i).ToString() + "*" + "P" + (i + _numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 1) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }
        public string dphi3ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q3controlpoints;
            for (int i = 0; i < valueVector.Count - 1; i++)
            {
                str += (_numControlPoints - 1 - i).ToString() + "*" + "P" + (i + 2 * _numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 1) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }

        public string ddphi1ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q1controlpoints;
            for (int i = 0; i < valueVector.Count - 2; i++)
            {
                str += ((_numControlPoints - 1 - i) * (_numControlPoints - 1 - i - 1)).ToString() + "*" + "P" + i.ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 2) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }
        public string ddphi2ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q2controlpoints;
            for (int i = 0; i < valueVector.Count - 2; i++)
            {
                str += ((_numControlPoints - 1 - i) * (_numControlPoints - 1 - i - 1)).ToString() + "*" + "P" + (i + _numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 2) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }
        public string ddphi3ToString()
        {
            string str = "";
            Vector<double> valueVector = _bezierMatrix * _q3controlpoints;
            for (int i = 0; i < valueVector.Count - 2; i++)
            {
                str += ((_numControlPoints - 1 - i) * (_numControlPoints - 1 - i - 1)).ToString() + "*" + "P" + (i + 2 * _numControlPoints).ToString() + "*" + "(theta^" + (_numControlPoints - 1 - i - 2) + ")" + "+";
            }
            str += "0";
            str = str.Replace("+-", "-");
            return str;
        }
    }
}
