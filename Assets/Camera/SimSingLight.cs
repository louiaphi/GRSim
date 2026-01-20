using JetBrains.Annotations;
using System.Text;
using Unity.VisualScripting;
using UnityEditor;
using UnityEngine;
using UnityEngine.XR;
using static UnityEditor.PlayerSettings;

public class SimSingLight : MonoBehaviour
{
    Matrix4x4 MetricTensor;
    Matrix4x4[] dMetricTensor;
    Matrix4x4[] Christoffelsymbols;
    Matrix4x4 diag;

    float M = 1;
    float a = 0.1f;

    void Start()
    {
        diag = Matrix4x4.identity;
        diag[0, 0] = -1;
        MetricTensor = Matrix4x4.identity;
        dMetricTensor = new Matrix4x4[4];
        Christoffelsymbols = new Matrix4x4[4];

        //Test();
        float r = CalcR(new Vector4(0, 1, 0, 0), a);
        CalcMetricTensor(new Vector4(0, 1, 0, 0), r, a, M);
        Debug.Log(r);
        LogTensor(MetricTensor);
        CopyMatrix(MetricTensor);
    }

    public void Test()
    {
        LightRay testLight1 = new LightRay(new Vector4(0, 100, 100, 0), new Vector4(1, 1, 0, 0));
        LightRay testLight2 = StepLight(testLight1, 1);
        LogRay(testLight2);
    }

    public void LogRay(LightRay ray)
    {
        Debug.Log("pt:");
        Debug.Log(ray.pos.x);
        Debug.Log("px:");
        Debug.Log(ray.pos.y);
        Debug.Log("py:");
        Debug.Log(ray.pos.z);
        Debug.Log("pz:");
        Debug.Log(ray.pos.w);
        Debug.Log("kt");
        Debug.Log(ray.k.x);
        Debug.Log("kx");
        Debug.Log(ray.k.y);
        Debug.Log("ky");
        Debug.Log(ray.k.z);
        Debug.Log("kz");
        Debug.Log(ray.k.w);
    }

    public void LogTensor(Matrix4x4 m)
    {
        Debug.Log("(" + m[0, 0] + ", " + m[0, 1] + ", " + m[0, 2] + ", " + m[0, 3] + ")");
        Debug.Log("(" + m[1, 0] + ", " + m[1, 1] + ", " + m[1, 2] + ", " + m[1, 3] + ")");
        Debug.Log("(" + m[2, 0] + ", " + m[2, 1] + ", " + m[2, 2] + ", " + m[2, 3] + ")");
        Debug.Log("(" + m[3, 0] + ", " + m[3, 1] + ", " + m[3, 2] + ", " + m[3, 3] + ")");
    }


    public static void CopyMatrix(Matrix4x4 m) //function written by ChatGPT for debugging
    {
        var sb = new StringBuilder();
        for (int row = 0; row < 4; row++)
        {
            sb.AppendLine(
                $"{m[row, 0]}, {m[row, 1]}, {m[row, 2]}, {m[row, 3]}"
            );
        }

        GUIUtility.systemCopyBuffer = sb.ToString();
        Debug.Log("Matrix copied to clipboard");
    }

    public float CalcR(Vector4 pos, float a)
    {
        float r = Mathf.Sqrt(CalcR2(pos, a));
        return r;
    }

    public struct LightRay
    {
        public LightRay(Vector4 pos, Vector4 k)
        {
            this.pos = pos;
            this.k = k;
        }
        public Vector4 pos;
        public Vector4 k; // Richtung
        
    }

    public LightRay StepLight (LightRay Y, float h)
    {
        //Beschleunigung
        //k_n wird als LightRay gespeichert, da es die selbe Form hat, ist aber anderer Natur (x->k,k->a)
        LightRay k1 = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        LightRay k2 = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        LightRay k3 = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        LightRay k4 = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        //k1
        k1 = CalculateA(Y);
        //k2
        LightRay k2TrialState = new LightRay(); //make new approximation
        for (int mu = 0; mu < 4; mu++)
        {
            k2TrialState.pos[mu] = Y.pos[mu] + (h / 2) * k1.pos[mu]; //k1.pos actually k
            k2TrialState.k[mu] = Y.k[mu] + (h/2) * k1.k[mu]; //k1.k actually a
        }
        k2 = CalculateA(k2TrialState);
        //k3
        LightRay k3TrialState = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        for (int mu = 0; mu < 4; mu++)
        {
            k3TrialState.pos[mu] = Y.pos[mu] + (h / 2) * k2.pos[mu]; //k2.pos actually k
            k3TrialState.k[mu] = Y.k[mu] + (h / 2) * k2.k[mu]; //k2.k actually a
        }
        k3 = CalculateA(k3TrialState);
        //k4
        LightRay k4TrialState = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        for (int mu = 0; mu < 4; mu++)
        {
            k4TrialState.pos[mu] = Y.pos[mu] + h * k3.pos[mu]; //k2.pos actually k
            k4TrialState.k[mu] = Y.k[mu] + h * k3.k[mu]; //k2.k actually a
        }
        k4 = CalculateA(k4TrialState);
        float posX = Y.pos.x + (h / 6) * (k1.pos.x + 2 * k2.pos.x + 2 * k3.pos.x + k4.pos.x); //Ich bereue es leicht Copilot deaktiviert zu haben......
        float posY = Y.pos.y + (h / 6) * (k1.pos.y + 2 * k2.pos.y + 2 * k3.pos.y + k4.pos.y);
        float posZ = Y.pos.z + (h / 6) * (k1.pos.z + 2 * k2.pos.z + 2 * k3.pos.z + k4.pos.z);
        float posW = Y.pos.w + (h / 6) * (k1.pos.w + 2 * k2.pos.w + 2 * k3.pos.w + k4.pos.w);
        float kX = Y.k.x + (h / 6) * (k1.k.x + 2 * k2.k.x + 2 * k3.k.x + k4.k.x);
        float kY = Y.k.y + (h / 6) * (k1.k.y + 2 * k2.k.y + 2 * k3.k.y + k4.k.y);
        float kZ = Y.k.z + (h / 6) * (k1.k.z + 2 * k2.k.z + 2 * k3.k.z + k4.k.z);
        float kW = Y.k.w + (h / 6) * (k1.k.w + 2 * k2.k.w + 2 * k3.k.w + k4.k.w);
        LightRay Y_new = new LightRay(new Vector4(posX, posY, posZ, posW), new Vector4(kX, kY, kZ, kW));
        Y_new.k.x = AdjustForConstraints(Y_new);
        return Y_new;
    }

    public LightRay CalculateA (LightRay aprxLight) //position, k -> k_n(k, a)
    {
        CalculateChristoffelsymbols(aprxLight.pos, Mathf.Sqrt(CalcR2(aprxLight.pos, a)), a, M);
        LightRay Output = new LightRay(new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
        for (int mu = 0; mu < 4; mu++)
        {
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    Output.k[mu] += -1 * Christoffelsymbols[mu][alpha, beta] * aprxLight.k[alpha] * aprxLight.k[beta];
                }
            }
        }
        Output.pos = aprxLight.k;
        return Output;
    }

    public  float AdjustForConstraints(LightRay Light) //warum bin ich so schrecklich im Namen geben???
    {
        CalcMetricTensor(Light.pos, Mathf.Sqrt(CalcR2(Light.pos, a)), a, M);
        float pd2 = 0; //p /2
        for (int i = 0; i < 4; i++)
        {
            pd2 += MetricTensor[0, i] * Light.k[i]; 
        }
        float Q = 0;
        for (int i = 0;i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                Q += MetricTensor[i, j] * Light.k[i] * Light.k[j];
            }
        }
        Q *= MetricTensor[0, 0];
        float k0 = (pd2 + Mathf.Sqrt(pd2 * pd2 - Q)) / (MetricTensor[0, 0]);
        return k0 * k0;
    }

    #region Metric Tensor Calculation

    public void CalcDelMetricTensor (Vector4 position, float r, float a, float M)
    {
        Vector4 L = CalcL(position, r, a);
        Vector4 delR = CalcDelR(position, a);
        Vector4[] delL = CalcDelL(position, r, delR, a);
        Vector4 delH = CalcDelH(position, r, M, a);
        float H = CalcH(r, position.w, M, a); //actually z, but (t, x, y, z) -> (x, y, z, w)

        for (int lambda = 1; lambda < 4; lambda++)
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
    public void CalcMetricTensor (Vector4 position, float r, float a, float M)
    {
        float z = position.w; //ahhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
        float H = CalcH(r, z, M, a);
        Vector4 l = CalcL(position, r, a);
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                MetricTensor[mu, nu] = diag[mu, nu] + 2 * H * l[mu] * l[nu];
            }
        }
    }

    public void CalculateChristoffelsymbols(Vector4 pos, float r, float a, float M)
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
    Vector4 CalcL(Vector4 pos, float r, float a)
    {
        Vector4 L = new Vector4(-1, (r * pos.y + a * pos.z) / (r * r + a * a), (r * pos.z - a * pos.y) / (r * r + a * a), pos.w / r);
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

    Vector4 CalcDelH(Vector4 pos, float r, float M, float a)
    {
        float delH_R = CalcDelH_R(r, pos.w, M, a);
        float delH_Z = CalcDelH_Z(r, pos.w, M, a);
        Vector4 dr = CalcDelR(pos, a);
        float delH_X = delH_R * dr.y;
        float delH_Y = delH_R * dr.z;
        float delH_Z_final = delH_R * dr.w + delH_Z;
        return new Vector4(0, delH_X, delH_Y, delH_Z_final);
    }

    #region Derivatives of L

    public Vector4[] CalcDelL(Vector4 pos, float r, Vector4 delR, float a) //
    {
        int[] KroX = new int[4] { 0, 1, 0, 0 }; //Kroenecker X
        int[] KroY = new int[4] { 0, 0, 1, 0 };
        int[] KroZ = new int[4] { 0, 0, 0, 1 };
        Vector4 delLT = new Vector4(0, 0, 0, 0); //partial derivative in all dircetions of the time component of the null covector is zero
        Vector4 delLX = new Vector4(0, 0, 0, 0); //just initialization
        Vector4 delLY = new Vector4(0, 0, 0, 0);
        Vector4 delLZ = new Vector4(0, 0, 0, 0);

        float FactorX = (pos.y * (r * r + a * a) - 2 * r * (r * pos.y + a * pos.z)) / (Mathf.Pow(r * r + a * a, 2));
        for (int dim = 1; dim < 4; dim++)
        {
            delLX[dim] = (r * KroX[dim] + a * KroY[dim]) / (r * r + a * a) + FactorX * delR[dim]; 
        }
        float FactorY = (pos.z * (r * r + a * a) - 2 * r * (r * pos.z - a * pos.y)) / (Mathf.Pow(r * r + a * a, 2));
        for (int dim = 1; dim < 4; dim++)
        {
            delLY[dim] = (r * KroY[dim] - a* KroX[dim]) / (r*r + a*a) + FactorY * delR[dim];
        }
        for (int dim = 1; dim < 4; dim++)
        {
            delLZ[dim] = KroZ[dim] - pos.w / (r*r) * delR[dim];
        }
        Vector4[] Output = new Vector4[4] { delLT, delLX, delLY, delLZ };
        return Output;

    }

    #endregion

    public float CalcR2(Vector4 pos, float a)
    {
        float rho2 = pos.y + pos.z + pos.w;
        return 0.5f * (rho2 - Mathf.Pow(a, 2) + Mathf.Sqrt(Mathf.Pow(rho2 - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.w, 2)));
    }

    public Vector4 CalcDelR2(Vector4 pos, float a)
    {
        float a2 = Mathf.Pow(a, 2);
        float rho2 = pos.y * pos.y + pos.z * pos.z + pos.w * pos.w;
        float rho = Mathf.Sqrt(rho2);
        float delRho2 = 2 * pos.y + 2 * pos.z + 2 * pos.w;
        float Delta = Mathf.Pow(pos.y + pos.z + pos.w - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.w, 2);
        float delDeltaX = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.y;
        float delDeltaY = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.z;
        float delDeltaZ = 2 * (rho2 - Mathf.Pow(a, 2)) * pos.w + 8 * Mathf.Pow(a, 2) * pos.w;
        float sqrtDelta = Mathf.Sqrt(Delta);
        float delR2X = 0.5f * (2 * pos.y + (2 * (rho2 - a2) * pos.y) / (sqrtDelta));
        float delR2Y = 0.5f * (2 * pos.z + (2 * (rho2 - a2) * pos.z) / (sqrtDelta));
        float delR2Z = 0.5f * (2 * pos.w + (2 * (rho2 - a2) * pos.w + 4 * a2 * pos.w) / (sqrtDelta));
        return new Vector4(0, delR2X, delR2Y, delR2Z);
    }

    public Vector4 CalcDelR(Vector4 pos, float a)
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
