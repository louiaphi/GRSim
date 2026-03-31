using Unity.VisualScripting;
using UnityEngine;

public class MainRender : MonoBehaviour
{

    Matrix4x4 MetricTensor;
    Matrix4x4 diag;

    public Vector4 CameraPosition = new Vector4(0, 100, 0, 0);
    public Vector3 CameraRotation;
    public float FOV;
    public Vector2Int MonitorSize = new Vector2Int(1080, 1920);

    Matrix4x4 MetricTensorAtCam;
    Matrix4x4 localTetradAtCam;

    public float M = 1;
    public float a = 0.1f;

    public float StepSize = 0.2f;
    public float preciseSteps = 0.1f;

    public float Rsoi;

    public int maxSteps;
    public float margin = 1;
    public float SuckInMargin;

    public ComputeShader fieldCS;
    RenderTexture target;

    void Start()
    {
        diag = Matrix4x4.identity;
        diag[0, 0] = -1;
        MetricTensor = Matrix4x4.identity;
    }

    void Dispatch()
    {
        if (target == null)
        {
            target = new RenderTexture(MonitorSize.x, MonitorSize.y, 0);
            target.enableRandomWrite = true;
            target.filterMode = FilterMode.Point;
            target.Create();
        }

        int kernel = fieldCS.FindKernel("CSMain");

        fieldCS.SetTexture(kernel, "Result", target);
        fieldCS.SetVector("CameraPosition", CameraPosition);
        fieldCS.SetVector("CameraRotation", CameraRotation);
        fieldCS.SetFloat("FOV", FOV);
        fieldCS.SetInts("MonitorSize", MonitorSize.x, MonitorSize.y);
        fieldCS.SetMatrix("MetricTensorAtCam", MetricTensorAtCam);
        fieldCS.SetMatrix("localTetradAtCam", localTetradAtCam);
        fieldCS.SetFloat("M", M);
        fieldCS.SetFloat("a", a);
        fieldCS.SetFloat("StepSize", StepSize);
        fieldCS.SetFloat("preciseSteps", preciseSteps);
        fieldCS.SetFloat("preciseSteps", preciseSteps);
        fieldCS.SetFloat("margin", margin);
        fieldCS.SetInt("maxSteps", maxSteps);
        fieldCS.SetFloat("Rsoi", Rsoi);
        fieldCS.SetFloat("SuckInMargin", SuckInMargin);
        fieldCS.Dispatch(
            kernel,
            Mathf.CeilToInt(MonitorSize.x / 8f),
            Mathf.CeilToInt(MonitorSize.y / 8f),
            1
        );
        GetComponent<Renderer>().material.mainTexture = target;
    }

    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.P))
        {
            Debug.Log("Started");
            CalcMetricTensor(CameraPosition, CalcR(CameraPosition, a), a, M);
            MetricTensorAtCam = MetricTensor;
            localTetradAtCam = localTetrad(CameraPosition);
            localTetradAtCam = rotateLocalTetrad(localTetradAtCam, CameraRotation);
            Dispatch();
            Debug.Log("finished");
        }
    }
    public float CalcR2(Vector4 pos, float a)
    {
        float rho2 = pos.y * pos.y + pos.z * pos.z + pos.w * pos.w;
        return 0.5f * (rho2 - Mathf.Pow(a, 2) + Mathf.Sqrt(Mathf.Pow(rho2 - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.w, 2)));
    }

    public float CalcR(Vector4 pos, float a)
    {
        float r = Mathf.Sqrt(CalcR2(pos, a));
        return r;
    }

    public void CalcMetricTensor(Vector4 position, float r, float a, float M)
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

    public Matrix4x4 localTetrad(Vector4 pos)
    {
        Vector4 u = new Vector4(1 / Mathf.Sqrt(-MetricTensorAtCam[0, 0]), 0, 0, 0);
        Vector4 v1 = new Vector4(0, 1, 0, 0);
        Vector4 v2 = new Vector4(0, 0, 1, 0);
        Vector4 v3 = new Vector4(0, 0, 0, 1);
        Matrix4x4 v = new Matrix4x4(u, v1, v2, v3);
        Matrix4x4 sqv = Matrix4x4.zero; //squiggely v~
        for (int BaseVector = 1; BaseVector < 4; BaseVector++)
        {
            for (int Component = 0; Component < 4; Component++)
            {
                float Q = 0;
                float D = 0;
                for (int alpha = 0; alpha < 4; alpha++)
                {
                    for (int beta = 0; beta < 4; beta++)
                    {
                        Q += MetricTensorAtCam[alpha, beta] * v[alpha, BaseVector] * u[beta];
                        D += MetricTensorAtCam[alpha, beta] * u[alpha] * u[beta];
                    }
                }
                sqv[Component, BaseVector] = v[Component, BaseVector] - Q / D * u[Component]; //check if indeces are the right way around
            }
        }
        Vector4 e0 = u;
        Vector4 e1 = new Vector4();
        Vector4 e2 = new Vector4();
        Vector4 e3 = new Vector4();
        Vector4 sqvs2 = new Vector4(); //squiggely v strich: ~v'_2
        Vector4 sqvs3 = new Vector4(); //squiggely v strich ~v'_3
        for (int dim = 0; dim < 4; dim++) //calc e1
        {
            float sum = 0;
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    sum += MetricTensorAtCam[alpha, beta] * sqv[alpha, 1] * sqv[beta, 1];
                }
            }
            e1[dim] = sqv[dim, 1] / Mathf.Sqrt(sum);
        }
        for (int dim = 0; dim < 4; dim++) //calc sqvs2
        {
            float sum = 0;
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    sum += MetricTensorAtCam[alpha, beta] * sqv[alpha, 2] * e1[beta];
                }
            }
            sqvs2[dim] = sqv[dim, 2] - sum * e1[dim];
        }
        for (int dim = 0; dim < 4; dim++)
        {
            float sum = 0;
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    sum += MetricTensorAtCam[alpha, beta] * sqvs2[alpha] * sqvs2[beta];
                }
            }
            e2[dim] = sqvs2[dim] / Mathf.Sqrt(sum);
        }
        for (int dim = 0; dim < 4; dim++)
        {
            float sum1 = 0;
            float sum2 = 0;
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    sum1 += MetricTensorAtCam[alpha, beta] * sqv[alpha, 3] * e1[beta];
                    sum2 += MetricTensorAtCam[alpha, beta] * sqv[alpha, 3] * e2[beta];
                }
            }
            sqvs3[dim] = sqv[dim, 3] - sum1 * e1[dim] - sum2 * e2[dim];
        }
        for (int dim = 0; dim < 4; dim++)
        {
            float sum = 0;
            for (int alpha = 0; alpha < 4; alpha++)
            {
                for (int beta = 0; beta < 4; beta++)
                {
                    sum += MetricTensorAtCam[alpha, beta] * sqvs3[alpha] * sqvs3[beta];
                }
            }
            e3[dim] = sqvs3[dim] / Mathf.Sqrt(sum);
        }
        Matrix4x4 e = new Matrix4x4(e0, e1, e2, e3);
        return e;

    }

    public Matrix4x4 rotateLocalTetrad(Matrix4x4 Tetrad, Vector3 Angles)
    {
        Matrix4x4 Rtest = Matrix4x4.Rotate(Quaternion.Euler(Angles));  // :skull: ... why did i listen to ChatGPT and deactivate my brain and write down all the rotation matrices myself
        Matrix4x4 R = Matrix4x4.zero;
        for (int i = 1; i < 4; i++)
        {
            for (int j = 1; j < 4; j++)
            {
                R[i, j] = Rtest[i - 1, j - 1];
            }
        }
        R[0, 0] = 1; //do not disturb time
        Matrix4x4 e = new Matrix4x4();
        for (int mu = 0; mu < 4; mu++)
        {
            for (int a = 0; a < 4; a++)
            {
                float sum = 0;
                for (int b = 0; b < 4; b++)
                {
                    sum += Tetrad[mu, b] * R[b, a];
                }
                e[mu, a] = sum;
            }
        }
        return e;
    }

}

