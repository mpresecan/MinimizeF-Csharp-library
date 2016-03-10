using System;

namespace MinimizeFCom
{

    public class Vector : Matrix
    {
        public Vector(int size) : base(size, 1) { }
        public Vector(string path) : base(path) { }
        public Vector(double[] val) : base(val) { }

        public int Length { get { return this.Nrow; } }

        public double this[int element]
        {
            get
            {
                if (this.matrix.GetLength(1) == 1)  //normal matrix nx1
                    return this.matrix[element - 1, 0];
                else                                //transposed matrix 1xn
                    return this.matrix[0, element - 1];
            }
            set
            {
                if (this.matrix.GetLength(1) == 1)
                    this.matrix[element - 1, 0] = value;
                else
                    this.matrix[0, element - 1] = value;
            }
        }

        /// <summary>
        /// Make new instance of Vector
        /// </summary>
        /// <returns></returns>
        //public new Vector Clone()
        //{
        //    Vector ret = new Vector(this.Length);
        //    for (int i = 0; i < this.Length; i++)
        //    {
        //        ret[i + 1] = this.matrix[i, 0];
        //    }
        //    return ret;
        //}


        // Overload operator '*' for Vector
        public static Vector operator *(double a, Vector v)
        {
            Vector ret = new Vector(v.Length);
            for (int i = 1; i <= v.Length; i++)
            {
                ret[i] = a * v[i];
            }

            return ret;
        }

        // Overload operator '/' for Vector
        public static Vector operator /(Vector v, double a)
        {
            if (a == 0) throw new DivideByZeroException();
            Vector ret = new Vector(v.Length);
            for (int i = 1; i <= v.Length; i++)
            {
                ret[i] = v[i] / a;
            }

            return ret;
        }

        // Overload operator '-' for Vector
        public static Vector operator -(Vector v1, Vector v2)
        {
            if (!Matrix.IsSameDimensions(v1, v2))
                throw new Exception("Matrixes are not the same dimmensions");

            Vector ret = new Vector(v1.Length);
            for (int i = 1; i <= v1.Length; i++)
            {
                ret[i] = v1[i] - v2[i];
            }

            return ret;
        }

        // Overload operator '+' for Vector
        public static Vector operator +(Vector v1, Vector v2)
        {
            if (!Matrix.IsSameDimensions(v1, v2))
                throw new Exception("Matrixes are not the same dimmensions");

            Vector ret = new Vector(v1.Length);
            for (int i = 1; i <= v1.Length; i++)
            {
                ret[i] = v1[i] + v2[i];
            }

            return ret;
        }

        /// <summary>
        /// String representation of Vector in '[',']'
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string ret = "[ ";
            for (int i = 0; i < this.matrix.GetLength(0) - 1; i++)
            {
                ret += this.matrix[i, 0] + ", ";
            }

            ret += this.matrix[this.matrix.GetLength(0) - 1, 0] + "]";
            return ret;
        }

        /// <summary>
        /// Calculate the norm of the Vector
        /// </summary>
        /// <returns></returns>
        public double Norm()
        {
            double sum = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                sum += Math.Pow(matrix[i, 0], 2);
            }
            return Math.Pow(sum, 0.5);
        }

        /// <summary>
        /// String representation of Vector without '[',']'
        /// </summary>
        /// <returns></returns>
        public string ToOutputString()
        {
            string output = "";
            for (int i = 1; i <= this.Length; i++)
            {
                output += this[i] + " ";
            }

            return output;
        }
    }
}
