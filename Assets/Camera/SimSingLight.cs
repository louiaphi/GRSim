using JetBrains.Annotations;
using UnityEditor;
using UnityEngine;
using UnityEngine.XR;

public class SimSingLight : MonoBehaviour
{
    Matrix4x4 MetricTensor;
    Matrix4x4[] dMetricTensor;
    Matrix4x4[] Christoffelsymbols;
    Matrix4x4 diag;

    float M = 1.0f;
    float a = 0.5f;

    void Start()
    {
        diag = Matrix4x4.identity;
        diag[0, 0] = -1;
        dMetricTensor = new Matrix4x4[4];
        Christoffelsymbols = new Matrix4x4[4];
    }

    #region Metric Tensor Calculation

    public Vector4 RK4(Vector3 pos, float StepSize) //4D space-Time pos not needed because integrating with affine partameter lambda
    {
        return new Vector4 ( 0, 0, 0, 0);
    }

    public void CalcDelMetricTensor(Vector3 pos, float r, float a, float M)
    {
        //Matrix4x4 delMetricTensorDim1 = MatrixScalarMul(OuterProduct(CalcL(pos, r, a), CalcL(pos, r, a)), 2 * CalcDelH(pos, r, M, a)[dim - 1]); //-1 weil dim 0 = Zeit //alt egal
        Vector4 L = CalcL(pos, r, a);
        Vector4 delR = CalcDelR(pos, a);
        Vector4[] delL = CalcDelL(pos, r, delR, a);
        Vector4 delH = CalcDelH(pos, r, M, a);
        float H = CalcH(r, pos.z, M, a);

        for (int lambda = 1; lambda  < 4; lambda++)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                for (int nu = 0; nu < 4; nu++)
                {
                    dMetricTensor[lambda][mu, nu] = 2 * (delH[lambda] * L[mu] * L[nu] + H * delL[lambda][mu] * L[nu] + H * delL[lambda][nu] * L[mu]);
                }
            }
        }
    }
    public void CalcMetricTensor(Vector3 position, float r, float a, float M)
    {
        float x = position.x;
        float y = position.y;
        float z = position.z;
        float x2 = Mathf.Pow(x, 2);
        float y2 = Mathf.Pow(y, 2);
        float z2 = Mathf.Pow(z, 2);
        float H = CalcH(r, z, M, z);
        float r2 = Mathf.Pow(r, 2);
        MetricTensor = new Matrix4x4(
            new Vector4(1, x / r, y / r, z / r),
            new Vector4(x / r, x2 / r2, x * y / r2, x * z / r2),
            new Vector4(y / r, x * y / r2, y2 / r2, y * z / r2),
            new Vector4(z / r, x * z / r2, y * z / r2, z2 / r2)
        );
        MetricTensor = MatrixScalarMul(MetricTensor, H);
        MetricTensor = MatrixAdd(MetricTensor, diag);
    }

    public void CalculateChristoffelsymbols(Vector3 pos, float r, float a, float M)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            CalcMetricTensor(pos, r, a, M);
            CalcDelMetricTensor(pos, r, M, a);
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    Christoffelsymbols[mu][alpha, beta] = 0;
                    for (int nu = 0; nu < 4; nu++)
                    {
                        Christoffelsymbols[mu][alpha, beta] += MetricTensor[mu, nu] * (dMetricTensor[alpha][nu, beta] + dMetricTensor[beta][nu, alpha] - dMetricTensor[nu][alpha, beta]);
                    }
                    Christoffelsymbols[mu][alpha, beta] = 0.5f * Christoffelsymbols[mu][alpha, beta];
                }
            }
        }
    }


    #region Auxiliary Functions for Metric Tensor
    Vector4 CalcL(Vector3 pos, float r, float a)
    {
        Vector4 L = new Vector4(1, (r * pos.x + a * pos.y) / (r * r + a * a), (r * pos.y - a * pos.x) / (r * r + a * a), pos.z / r);
        return L;
    }

    float CalcH(float r, float z, float M, float a) //static strength of deviation from flat Minkowski Spacetime space (H = 0 -> Minkowski) at Horizon H = infinty
    {
        float H = (M * Mathf.Pow(r, 3)) / (Mathf.Pow(r, 4) + Mathf.Pow(a, 2) * Mathf.Pow(z, 2));
        return H;
    }

    float CalcDelH_R(float r, float z, float M, float a)
    {
        float numerator = M * (3 * Mathf.Pow(r, 2) * (Mathf.Pow(r, 4) + Mathf.Pow(a, 2) * Mathf.Pow(z, 2)) - Mathf.Pow(r, 3) * (4 * Mathf.Pow(r, 3)));
        float denominator = Mathf.Pow(Mathf.Pow(r, 4) + Mathf.Pow(a, 2) * Mathf.Pow(z, 2), 2);
        return numerator / denominator;
    }

    float CalcDelH_Z(float r, float z, float M, float a)
    {
        float numerator = -2 * M * Mathf.Pow(r, 3) * Mathf.Pow(a, 2) * z;
        float denominator = Mathf.Pow(Mathf.Pow(r, 4) + Mathf.Pow(a, 2) * Mathf.Pow(z, 2), 2);
        return numerator / denominator;
    }

    Vector4 CalcDelH(Vector3 pos, float r, float M, float a)
    {
        float delH_R = CalcDelH_R(r, pos.z, M, a);
        float delH_Z = CalcDelH_Z(r, pos.z, M, a);
        Vector3 dr = CalcDelR(pos, a);
        float delH_X = delH_R * dr.x;
        float delH_Y = delH_R * dr.y;
        float delH_Z_final = delH_R * dr.z + delH_Z;
        return new Vector4(0, delH_X, delH_Y, delH_Z_final);
    }

    #region Derivatives of L

    public Vector4[] CalcDelL(Vector3 pos, float r, Vector4 delR, float a) //
    {
        int[] KroX = new int[4] { 0, 1, 0, 0 }; //Kroenecker X
        int[] KroY = new int[4] { 0, 0, 1, 0 };
        int[] KroZ = new int[4] { 0, 0, 0, 1 };
        Vector4 delLT = new Vector4(0, 0, 0, 0); //partial derivative in all dircetions of the time component of the null covector is zero
        Vector4 delLX = new Vector4(0, 0, 0, 0); //just initialization
        Vector4 delLY = new Vector4(0, 0, 0, 0);
        Vector4 delLZ = new Vector4(0, 0, 0, 0);

        float FactorX = (pos.x * (r * r + a * a) - 2 * r * (r * pos.x + a * pos.y)) / (Mathf.Pow(r * r + a * a, 2));
        for (int dim = 1; dim < 4; dim++)
        {
            delLX[dim] = (r * KroX[dim] + a * KroY[dim]) / (r * r + a * a) + FactorX * delR[dim]; 
        }
        float FactorY = (pos.y * (r * r + a * a) - 2 * r * (r * pos.y - a * pos.x)) / (Mathf.Pow(r * r + a * a, 2));
        for (int dim = 1; dim < 4; dim++)
        {
            delLY[dim] = (r * KroY[dim] - a* KroX[dim]) / (r*r + a*a) + FactorY * delR[dim];
        }
        for (int dim = 1; dim < 4; dim++)
        {
            delLZ[dim] = KroZ[dim] - pos.z / (r*r) * delR[dim];
        }
        Vector4[] Output = new Vector4[4] { delLT, delLX, delLY, delLZ };
        return Output;

    }

    //public Vector4 CalcDelL_x (Vector3 pos, float r, float M, float a) //derivatives of l_x in respect to x, y, z and r
    //{
    //   float DelLx_X = r / (Mathf.Pow(r, 2) + Mathf.Pow(a, 2));
    //   float DelLx_Y = a / (Mathf.Pow(r, 2) + Mathf.Pow(a, 2));
    //   float DelLx_Z = 0f;
    //   float DelLx_R = (pos.x * (Mathf.Pow(r, 2) + Mathf.Pow(a, 2)) - (r * pos.x + a * pos.y) * 2 * r) / Mathf.Pow((Mathf.Pow(r, 2) + Mathf.Pow(a, 2)), 2);
    //   return new Vector4(DelLx_X, DelLx_Y, DelLx_Z, DelLx_R);
    //}
    //
    //public Vector4 CalcDelL_y(Vector3 pos, float r, float M, float a) //derivatives of l_y in respect to x, y, z and r
    //{
    //    float DelLy_X = -a / (Mathf.Pow(r, 2) + Mathf.Pow(a, 2));
    //    float DelLy_Y = r / (Mathf.Pow(r, 2) + Mathf.Pow(a, 2));
    //    float DelLy_Z = 0f;
    //    float DelLy_R = (pos.y * (Mathf.Pow(r, 2) + Mathf.Pow(a, 2)) - (r * pos.y - a * pos.x) * 2 * r) / Mathf.Pow((Mathf.Pow(r, 2) + Mathf.Pow(a, 2)), 2);
    //    return new Vector4(DelLy_X, DelLy_Y, DelLy_Z, DelLy_R);
    //}
    //
    //public Vector4 CalcDelL_z(Vector3 pos, float r, float M, float a) //derivatives of l_z in respect to x, y, z and r
    //{
    //    float DelLz_X = 0f;
    //    float DelLz_Y = 0f;
    //    float DelLz_Z = 1 / r;
    //    float DelLz_R = -pos.z / (Mathf.Pow(r, 2));
    //    return new Vector4(DelLz_X, DelLz_Y, DelLz_Z, DelLz_R);
    //}

    #endregion

    public float CalcR2(Vector3 pos, float a)
    {
        float rho2 = pos.x + pos.y + pos.z;
        return 0.5f * (rho2 - Mathf.Pow(a, 2) + Mathf.Sqrt(Mathf.Pow(rho2 - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.z, 2)));
    }

    public Vector4 CalcDelR2(Vector3 pos, float a)
    {
        float a2 = Mathf.Pow(a, 2);
        float rho2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
        float rho = Mathf.Sqrt(rho2);
        float delRho2 = 2 * pos.x + 2 * pos.y + 2 * pos.z;
        float Delta = Mathf.Pow(pos.x + pos.y + pos.z - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.z, 2);
        float delDeltaX = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.x;
        float delDeltaY = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.y;
        float delDeltaZ = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.z + 8 * Mathf.Pow(a, 2) * pos.z;
        float sqrtDelta = Mathf.Sqrt(Delta);
        float delR2X = 0.5f * (2 * pos.x + (2 * (rho2 - a2) * pos.x) / (sqrtDelta));
        float delR2Y = 0.5f * (2 * pos.y + (2 * (rho2 - a2) * pos.y) / (sqrtDelta));
        float delR2Z = 0.5f * (2 * pos.z + (2 * (rho2 - a2) * pos.z + 4 * a2 * pos.z) / (sqrtDelta));
        return new Vector4(0, delR2X, delR2Y, delR2Z);
    }

    public Vector4 CalcDelR(Vector3 pos, float a)
    {
        Vector4 delR2 = CalcDelR2(pos, a);
        float R2 = CalcR2(pos, a);
        Vector4 delR = (1 / (2 * Mathf.Sqrt(R2)) * delR2);
        return delR;
    }

    #endregion

    #endregion


    public Matrix4x4 MatrixAdd(Matrix4x4 m1, Matrix4x4 m2)
    {
        Matrix4x4 result = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = m1[i, j] + m2[i, j];
            }
        }
        return result;
    }

    public Matrix4x4 MatrixScalarMul(Matrix4x4 m, float s)
    {
        Matrix4x4 result = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = m[i, j] * s;
            }
        }
        return result;
    }

    public Matrix4x4 OuterProduct(Vector4 v1, Vector4 v2)
    {
        Matrix4x4 result = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = v1[i] * v2[j];
            }
        }
        return result;
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
