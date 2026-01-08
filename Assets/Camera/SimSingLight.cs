using UnityEngine;

public class SimSingLight : MonoBehaviour
{
    Matrix4x4 MetricTensor;
    Matrix4x4 diag;

    float M = 1.0f;
    float a = 0.5f;

    void Start()
    {
        diag = Matrix4x4.identity;
        diag[0, 0] = -1;
    }

    public void CalcMetricTensor(Vector3 position)
    {
        float x = position.x;
        float y = position.y;
        float z = position.z;
        float x2 = Mathf.Pow(x, 2);
        float y2 = Mathf.Pow(y, 2);
        float z2 = Mathf.Pow(z, 2);
        float r = CalcR(position, a);
        float H = M / r;
        float r2 = Mathf.Pow(r, 2);
        MetricTensor = new Matrix4x4(
            new Vector4(1, x/r, y/r, z/r),
            new Vector4(x/r, x2/r2, x*y/r2, x*z/r2),
            new Vector4(y/r, x*y/r2, y2/r2, y*z/r2),
            new Vector4(z/r, x*z/r2, y*z/r2, z2/r2)
        );
        MetricTensor = MatrixScalarMul(MetricTensor, H);
        MetricTensor = MatrixAdd(MetricTensor, diag);
    }

    public float CalcR(Vector3 pos, float a)
    {
        float rho = pos.magnitude;
        return Mathf.Sqrt(0.5f * (Mathf.Pow(rho, 2) - Mathf.Pow(a, 2) + Mathf.Sqrt(Mathf.Pow(Mathf.Pow(rho, 2) - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.z, 2))));
    }
    
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
    // Update is called once per frame
    void Update()
    {
        
    }
}
