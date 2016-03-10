using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace MinimizeFCom
{
    public class Matrix
    {
        protected double[,] matrix;
        public int Nrow;
        public int Ncoll;
        public bool Squared { get; protected set; }
        public bool Transposed { get; protected set; } = false;

        //precision constant
        protected const double epsilon = 0.0000001;

        // ??
        enum All { };

        /// <summary>
        /// Create Matrix with size of 'row' X 'coll', all elements initialized on 0
        /// </summary>
        /// <param name="row">Numbers of rows</param>
        /// <param name="coll">Number of collumns</param>
        public Matrix(int row, int coll)
        {
            this.matrix = new double[row, coll];
            this.Nrow = row;
            this.Ncoll = coll;
            this.Squared = (row == coll) ? true : false;
        }

        /// <summary>
        /// Create squared Matrix size of 'size' X 'size'
        /// </summary>
        /// <param name="size">Size of squared matrix</param>
        public Matrix(int size)
        {
            this.matrix = new double[size, size];
            this.Nrow = size;
            this.Ncoll = size;
            this.Squared = true;
        }

        /// <summary>
        /// Create Matrix with values of one dimensional double array
        /// </summary>
        /// <param name="val">Array of values</param>
        public Matrix(double[] val)
        {
            this.matrix = new double[val.Length, 1];
            for (int i = 0; i < val.Length; i++)
            {
                this.matrix[i, 0] = val[i];
            }
            this.Nrow = val.Length;
            this.Ncoll = 1;
            this.Squared = false;
        }

        /// <summary>
        /// Create Matrix with values of two dimensional double array
        /// </summary>
        /// <param name="val">Two dimensional array of values</param>
        public Matrix(double[,] val)
        {
            this.matrix = new double[val.Length, 1];
            this.matrix = val;
            this.Nrow = val.GetLength(0);
            this.Ncoll = val.GetLength(1);
            this.Squared = (val.GetLength(0) == val.GetLength(1)) ? true : false;
        }

        /// <summary>
        /// Create Matrix from the file
        /// </summary>
        /// <param name="path">Path to the stored Matrix</param>
        public Matrix(string path)
        {
            List<string> lines = new List<string>();

            //download it from the path something like Matrix A("A.txt")
            using (StreamReader file = new StreamReader(path))
            {
                string line;
                while ((line = file.ReadLine()) != null)
                    lines.Add(line);
            }

            List<List<double>> matrixFromFile = new List<List<double>>();
            string[] numbers;
            foreach (var line in lines)
            {
                List<double> elements = new List<double>();
                numbers = line.Split(' ');
                foreach (var number in numbers)
                    elements.Add(Convert.ToDouble(number));
                matrixFromFile.Add(elements);
            }

            //create Matrix
            int row = matrixFromFile.Count();
            if (row != 0)
            {
                int coll = matrixFromFile.First().Count();
                Matrix temp = new Matrix(row, coll);

                this.matrix = new double[row, coll];

                for (int i = 0; i < row; i++)
                    for (int j = 0; j < coll; j++)
                        this.matrix[i, j] = matrixFromFile.ElementAt(i).ElementAt(j);

                this.Nrow = row;
                this.Ncoll = coll;
                this.Squared = (row == coll) ? true : false;
            }
            else
            {
                this.Nrow = 0;
                this.Ncoll = 0;
                this.Squared = false;
                this.matrix = null;
            }
        }

        /// <summary>
        /// Allocate a value in Matrix. One-base indexing.
        /// </summary>
        /// <param name="row"></param>
        /// <param name="coll"></param>
        /// <returns></returns>
        public double this[int row, int coll]
        {
            get
            {
                return this.matrix[row - 1, coll - 1];
            }
            set
            {
                this.matrix[row - 1, coll - 1] = value;
            }
        }

        /// <summary>
        /// Retreive specific row from the Matrix as Vector
        /// </summary>
        /// <param name="row">One-base index of a row</param>
        /// <returns>Row vector</returns>
        public Vector GetRowVector(int row)
        {
            Vector ret = new Vector(this.matrix.GetLength(1));
            for (int i = 0; i < this.matrix.GetLength(1); i++)
                ret[i + 1] = this.matrix[row - 1, i];

            return ret;
        }

        /// <summary>
        /// Retrieve specific collumn from Matrix as Vector
        /// </summary>
        /// <param name="coll">One-base indexe of a collumn vector</param>
        /// <returns></returns>
        public Vector GetCollumnVector(int coll)
        {
            Vector ret = new Vector(this.matrix.GetLength(0));
            for (int i = 0; i < this.matrix.GetLength(0); i++)
                ret[i + 1] = this.matrix[i, coll - 1];

            return ret;
        }

        /// <summary>
        /// Change entire row in Matrix by values in Vector
        /// </summary>
        /// <param name="rowN">One-base index of a row</param>
        /// <param name="rowVect">New values</param>
        public void SetRowVector(int rowN, Vector rowVect)
        {
            for (int i = 0; i < this.matrix.GetLength(1); i++)
                this.matrix[rowN - 1, i] = rowVect[i + 1];
        }

        /// <summary>
        /// Change entire collumn in Matrix by values in Vector
        /// </summary>
        /// <param name="collN">One-base index of a collumn</param>
        /// <param name="collVect">New values</param>
        public void SetCollumnVector(int collN, Vector collVect)
        {
            for (int i = 0; i < this.matrix.GetLength(0); i++)
                this.matrix[i, collN - 1] = collVect[i + 1];
        }

        // Overload operator '+'
        public static Matrix operator +(Matrix m1, Matrix m2)
        {
            if (!Matrix.IsSameDimensions(m1, m2))
                throw new Exception("Matrixes are not the same dimmensions");

            Matrix temp = new Matrix(m1.Nrow, m1.Ncoll);

            for (int i = 1; i <= m1.Nrow; i++)
            {
                for (int j = 1; j <= m1.Ncoll; j++)
                {
                    temp[i, j] = m1[i, j] + m2[i, j];
                }
            }
            return temp;
        }


        // Overload operator '-'
        public static Matrix operator -(Matrix m1, Matrix m2)
        {
            if (!Matrix.IsSameDimensions(m1, m2))
                throw new Exception("Matrixes are not the same dimmensions");

            Matrix temp = new Matrix(m1.Nrow, m1.Ncoll);

            for (int i = 1; i <= m1.Nrow; i++)
            {
                for (int j = 1; j <= m1.Ncoll; j++)
                {
                    temp[i, j] = m1[i, j] - m2[i, j];
                }
            }
            return temp;
        }

        // overloading operator '=='
        public static bool operator ==(Matrix m1, Matrix m2)
        {
            // Check reference equality
            if (object.ReferenceEquals(m1, m2)) return true;
            else if ((object.ReferenceEquals(m1, null) && !object.ReferenceEquals(m2, null)) || (!object.ReferenceEquals(m1, null) && object.ReferenceEquals(m2, null))) return false;

            if (!Matrix.IsSameDimensions(m1, m2))
                return false;

            // Check value equality
            for (int i = 1; i <= m1.Nrow; i++)
            {
                for (int j = 1; j <= m1.Ncoll; j++)
                {
                    if (Math.Abs(m1[i, j] - m2[i, j]) > Matrix.epsilon) return false;
                }
            }
            return true;
        }

        public static bool operator !=(Matrix m1, Matrix m2)
        {
            return !(m1 == m2);
        }

        // Overload operator '*'
        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            if (m1.Ncoll != m2.Nrow)
                throw new Exception("Matrixes do not fit the dimensions to multiply");

            Matrix temp = new Matrix(m1.Nrow, m2.Ncoll);

            for (int i = 1; i <= m1.Nrow; i++)
            {
                for (int j = 1; j <= m2.Ncoll; j++)
                {
                    temp[i, j] = 0;
                    for (int k = 1; k <= m2.Nrow; k++)
                    {
                        temp[i, j] += m1[i, k] * m2[k, j];
                    }
                }
            }
            return temp;
        }

        // Overload operator '*' with vector
        public static Vector operator *(Matrix m1, Vector b)
        {
            if (m1.Ncoll != b.Nrow)
                throw new Exception("Matrixes do not fit the dimensions to multiply");

            Matrix temp = new Matrix(m1.Nrow, b.Ncoll);

            for (int i = 1; i <= m1.Nrow; i++)
            {
                for (int j = 1; j <= b.Ncoll; j++)
                {
                    temp[i, j] = 0;
                    for (int k = 1; k <= b.Nrow; k++)
                    {
                        temp[i, j] += m1[i, k] * b[k, j];
                    }
                }
            }
            return temp.ConvertToVector();
        }

        // Overload operator '*' - matrix * scalar
        public static Matrix operator *(Matrix m, double scalar)
        {
            Matrix temp = new Matrix(m.Nrow, m.Ncoll);

            for (int i = 1; i <= m.Nrow; i++)
                for (int j = 1; j <= m.Ncoll; j++)
                    temp[i, j] = m[i, j] * scalar;

            return temp;
        }

        // Overload operator '*' - scalar * matrix
        public static Matrix operator *(double scalar, Matrix m)
        {

            Matrix temp = new Matrix(m.Nrow, m.Ncoll);

            for (int i = 1; i <= m.Nrow; i++)
                for (int j = 1; j <= m.Ncoll; j++)
                    temp[i, j] = m[i, j] * scalar;

            return temp;
        }

        /// <summary>
        /// Multiply each element in 'm1' with 'm2'
        /// </summary>
        /// <param name="m1">First multiplier</param>
        /// <param name="m2">Second multiplier</param>
        /// <returns></returns>
        public static Matrix SimpleMultply(Matrix m1, Matrix m2)
        {
            if (!Matrix.IsSameDimensions(m1, m2))
                throw new Exception("Matrixes are not the same dimmensions");

            Matrix temp = new Matrix(m1.Nrow, m2.Ncoll);

            for (int i = 1; i <= m1.Nrow; i++)
                for (int j = 1; j <= m1.Ncoll; j++)
                    temp[i, j] = m1[i, j] * m2[i, j];

            return temp;
        }

        /// <summary>
        /// Transpose Matrix
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public Matrix Transpose(Matrix m)
        {
            Matrix temp = new Matrix(m.Ncoll, m.Nrow);
            for (int i = 1; i <= m.Ncoll; i++)
                for (int j = 1; j <= m.Nrow; j++)
                    temp[i, j] = m[j, i];

            m.Transposed = m.Transposed ? false : true;
            return temp;
        }

        /// <summary>
        /// Check if two Matrixes are the same dimensions
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static bool IsSameDimensions(Matrix m1, Matrix m2)
        {
            if (m1.Nrow == m2.Nrow && m1.Ncoll == m2.Ncoll) return true;
            else return false;
        }

        /// <summary>
        /// Convert Matrix into Vector. 2D -> 1D
        /// </summary>
        /// <returns></returns>
        public Vector ConvertToVector()
        {
            Vector ret;
            if (this.matrix.GetLength(0) == 1)
            {
                ret = new Vector(matrix.GetLength(1));
                for (int i = 0; i < this.matrix.GetLength(1); i++)
                {
                    ret[i + 1] = this.matrix[0, i];
                }
            }
            else if (this.matrix.GetLength(1) == 1)
            {
                ret = new Vector(matrix.GetLength(0));
                for (int i = 0; i < this.matrix.GetLength(0); i++)
                {
                    ret[i + 1] = this.matrix[i, 0];
                }
            }
            else
            {
                ret = new Vector(matrix.GetLength(0) * matrix.GetLength(1));
                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                        ret[i + j] = matrix[i, j];
                    }
                }
            }
            return ret;
        }

        /// <summary>
        /// String representation of a Matrix
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string output = "";
            for (int i = 0; i < this.Nrow; i++)
            {
                for (int j = 0; j < this.Ncoll; j++)
                {
                    output += matrix[i, j].ToString();
                    if (j < this.Ncoll - 1) output += " ";
                }
                output += Environment.NewLine;
            }
            return output;
        }

        /// <summary>
        /// Clone a Matrix
        /// </summary>
        /// <returns></returns>
        public Matrix Clone()
        {
            //has to be done for every property in class
            Matrix temp = new Matrix(this.Nrow, this.Ncoll);
            for (int i = 1; i <= this.Nrow; i++)
                for (int j = 1; j <= this.Ncoll; j++)
                    temp[i, j] = this.matrix[i - 1, j - 1];

            temp.Nrow = this.Nrow;
            temp.Ncoll = this.Ncoll;
            if (this.Nrow == this.Ncoll) temp.Squared = true;
            else temp.Squared = false;

            return temp;
        }

        /// <summary>
        /// Performe substitution forward
        /// </summary>
        /// <param name="A">Set of linear equation matrix</param>
        /// <param name="b">Solution</param>
        /// <param name="order">Permutation matrix</param>
        /// <returns>New subscore Y</returns>
        public static Vector SubstitutionForward(Matrix A, Vector b, Matrix order = null)
        {
            if (!ReferenceEquals(null, order))
            {
                b = order * b;
                if (order.Nrow != A.Nrow || order.Ncoll != A.Ncoll) throw new Exception("Permutation matrix and matrix of supstitution are not the same");
            }

            if (!A.Squared) throw new Exception("Matrix is not squared");

            for (int i = 1; i <= b.Length - 1; i++) //pazi dali je ukljućeno ili ne <= ili <
            {
                for (int j = i + 1; j <= b.Length; j++) //također pazi
                {
                    b[j] -= A[j, i] * b[i];
                }
            }

            return b;
        }

        /// <summary>
        /// Performe substitution backward
        /// </summary>
        /// <param name="A">Set of linear equation matrix</param>
        /// <param name="b">Subscore Y</param>
        /// <param name="e">Precision</param>
        /// <returns>Solution Vector X</returns>
        public static Vector SubstitutionBackward(Matrix A, Vector b, double e = epsilon)
        {
            if (!A.Squared) throw new Exception("Matrix is not squared");

            for (int i = b.Length; i >= 1; i--) //također pazi
            {
                if (Math.Abs(A[i, i]) < e) throw new Exception("Dividing with number close to zero");
                b[i] /= A[i, i];
                for (int j = 1; j <= i - 1; j++)
                {
                    b[j] -= A[j, i] * b[i];
                }
            }
            return b;
        }

        /// <summary>
        /// Makes a LU decomposition if P atribute is set false. Makes LUP decomposition if P is set true. Default it makes LUP decomposition.
        /// Perform LUP decomposition. If 'P' is set to 'false', then it performs LU decomposition.
        /// </summary>
        /// <param name="A">Matrix to make a decomposition. It changes the matrix</param>
        /// <param name="P">If true (default) -> LUP || false -> LU</param>
        /// <param name="e">Precision work</param>
        /// <returns>Permutated Indentity Matrix for LUP decomposition</returns>
        public static Matrix LUP(Matrix A, bool P = true, double e = epsilon)
        {
            if (!A.Squared) throw new Exception("Matrix is not squared");

            Matrix Acopy = A.Clone();

            Matrix E = Matrix.Eye(Acopy.Nrow);

            for (int i = 1; i <= Acopy.Nrow - 1; i++)
            {
                if (P) //Pivoting enabled - aditionaly performing LUP decomposition
                {
                    int R = choosePivotElement(Acopy, i);
                    if (Math.Abs(Acopy[R, i]) < e) throw new Exception("Dividing with number close to zero");
                    if (R != i)
                    {
                        Vector buffer = Acopy.GetRowVector(R);
                        Vector oldBuff = Acopy.GetRowVector(i);
                        Acopy.SetRowVector(i, buffer);
                        Acopy.SetRowVector(R, oldBuff);

                        buffer = E.GetRowVector(R);
                        oldBuff = E.GetRowVector(i);
                        E.SetRowVector(i, buffer);
                        E.SetRowVector(R, oldBuff);
                    }
                }

                for (int j = i + 1; j <= Acopy.Nrow; j++)
                {
                    if (Math.Abs(Acopy[i, i]) < e) throw new Exception("Dividing with number close to zero");
                    Acopy[j, i] /= Acopy[i, i];
                    for (int k = i + 1; k <= Acopy.Nrow; k++)
                    {
                        Acopy[j, k] -= Acopy[j, i] * Acopy[i, k];
                    }
                }
            }
            return E;
        }

        private static int choosePivotElement(Matrix A, int i)
        {
            int theBiggest = i;
            for (int j = i + 1; j <= A.Nrow; j++)
                if (Math.Abs(A[j, i]) > Math.Abs(A[theBiggest, i])) theBiggest = j;
            return theBiggest;
        }

        /// <summary>
        /// Create Eye Matrix
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns></returns>
        public static Matrix Eye(int n)
        {
            Matrix E = new Matrix(n);
            for (int i = 1; i <= n; i++)
                E[i, i] = 1;
            return E;
        }

        /// <summary>
        /// Solve N equations with N variables
        /// </summary>
        /// <param name="A">Equetions matrix</param>
        /// <param name="b">Free memeber vector</param>
        /// <param name="LUP">Do LUP decomposition, if false do LU decomposition</param>
        /// <param name="e">Precision</param>
        /// <returns></returns>
        public static Vector SolveEquations(Matrix A, Vector b, bool LUP = true, double e = epsilon, bool print = false)
        {
            if (LUP && e == epsilon)
            {
                Matrix P = Matrix.LUP(A);
                b = Matrix.SubstitutionForward(A, b, P);
                b = Matrix.SubstitutionBackward(A, b);
            }
            else if (LUP && e != epsilon)
            {
                Matrix P = Matrix.LUP(A, true, e);
                b = Matrix.SubstitutionForward(A, b, P);
                b = Matrix.SubstitutionBackward(A, b, e);
            }
            else
            {
                Matrix.LUP(A, false);
                b = Matrix.SubstitutionForward(A, b);
                b = Matrix.SubstitutionBackward(A, b);
            }
            return b;
        }

        /// <summary>
        /// Inverse Matrix
        /// </summary>
        /// <param name="A"></param>
        /// <returns>Inverse </returns>
        public static Matrix Inverse(Matrix A)
        {
            if (A.Nrow != A.Ncoll) throw new Exception("Matrix is not squared");
            Matrix E = Matrix.Eye(A.Ncoll);
            Matrix inverse = new Matrix(A.Nrow);

            Matrix LUP = A.Clone();
            Matrix P = Matrix.LUP(LUP);
            for (int i = 1; i <= A.Ncoll; i++)
            {
                Vector b = E.GetCollumnVector(i);
                b = Matrix.SubstitutionForward(LUP, b, P);
                b = Matrix.SubstitutionBackward(LUP, b);
                inverse.SetCollumnVector(i, b);
            }

            return inverse;

        }

        /// <summary>
        /// Save Matrix to file
        /// </summary>
        /// <param name="path">Path with name and type of file</param>
        public void WriteToFile(string path)
        {
            using (StreamWriter sw = new StreamWriter(path))
            {
                sw.Write(this.ToString());
            }
        }
    }
}
