using Unity.VisualScripting;
using UnityEngine;

public class ScalarFieldVisualizer : MonoBehaviour
{
    public int width = 500;
    public int height = 500;

    public float scale = 1.0f;

    public bool UpdateByTime = true;
    public float UpdateTime;

    public float a;
    public float M;
    public float z;

    public float rim;

    public bool enableEventHorizon;
    public bool enableErgosphere;

    Texture2D texture;

    public float delT = 0;

    public float min = 0;
    public float max = 0;
    [Range(0.0f, 10.0f)]
    public float xgrob; //two sliders for better precision control
    [Range(0.0f, 1.0f)]
    public float xfine;
    public float xpoint;

    void Start()
    {

    }

    //float CalcH(float r, float z, float M, float a) //static strength of deviation from flat Minkowski Spacetime space (H = 0 -> Minkowski) at Horizon H = infinty
    //{
    //    float H = (M * Mathf.Pow(r, 3)) / (Mathf.Pow(r, 4) + Mathf.Pow(a, 2) * Mathf.Pow(z, 2));
    //    return H;
    //}
    //
    //public float CalcR2(Vector4 pos, float a)
    //{
    //    float rho2 = pos.y * pos.y + pos.z * pos.z + pos.w * pos.w;
    //    return 0.5f * (rho2 - Mathf.Pow(a, 2) + Mathf.Sqrt(Mathf.Pow(rho2 - Mathf.Pow(a, 2), 2) + 4 * Mathf.Pow(a, 2) * Mathf.Pow(pos.w, 2)));
    //}
    //
    //void GenerateTexture() // made by ChatGPT for debugging
    //{
    //    Color[] pixels = new Color[width * height];
    //
    //
    //    float[,] values = new float[width, height];
    //    for (int y = 0; y < height; y++)
    //    {
    //        for (int x = 0; x < width; x++)
    //        {
    //            float v = SampleField(x, y);
    //            values[x, y] = v;
    //            //min = Mathf.Min(min, v);
    //            //max = Mathf.Max(max, v);
    //        }
    //    }
    //
    //    float range = max - min;
    //
    //    for (int y = 0; y < height; y++)
    //    {
    //        for (int x = 0; x < width; x++)
    //        {
    //            if (values[x, y] == -1)
    //            {
    //                pixels[y * width + x] = new Color(1, 0, 0, 1f);
    //            }
    //            else if (values[x, y] == -2)
    //            {
    //                pixels[y * width + x] = new Color(0, 1, 0, 1f);
    //            }
    //            float normalized = (values[x, y] - min) / range;
    //            pixels[y * width + x] = new Color(normalized, normalized, normalized, 1f);
    //        }
    //    }
    //
    //    texture.SetPixels(pixels);
    //    texture.Apply();
    //}
    //
    //float SampleField(int x, int y)
    //{
    //    float r = Mathf.Sqrt(CalcR2(new Vector4(0, x - width/2, y-width/2, z), a));
    //    if (r < 0.1f && r > -0.1f)
    //    {
    //        return -1;
    //    }
    //    if (r < -0.1f)
    //    {
    //        return -2;
    //}
    //    r *= scale;
    //    float H = CalcH(r, z, M, a);
    //    Debug.Log("R = " + Mathf.Sqrt(CalcR2(new Vector4(0, xpoint, 0, 0), a)));
    //    Debug.Log("H = " + CalcH(CalcR2(new Vector4(0, xpoint, 0, 0), a), 0, M, a));
    //    return H;
    //}

    public ComputeShader fieldCS;
    RenderTexture target;

    public float CalcEHRadius(float M, float a)
    {
        float r = M + Mathf.Sqrt(M * M - a * a);
        return r;
    }


    public float RAtZ (float r, float z)
    {
        return (r * r) * (1 - (z * z) / (r * r));
    }


    void Dispatch() //made by ChatGPT for Fun yay (same as old code but on GPU as Compute Shader so i can get nicer pictures that are absolutly mandatory for good debugging)
    {
        if (target == null)
        {
            target = new RenderTexture(width, height, 0);
            target.enableRandomWrite = true;
            target.filterMode = FilterMode.Point;
            target.Create();
        }

        int kernel = fieldCS.FindKernel("CSMain");

        fieldCS.SetTexture(kernel, "Result", target);
        fieldCS.SetInt("width", width);
        fieldCS.SetInt("height", height);
        fieldCS.SetFloat("scale", scale);
        fieldCS.SetFloat("a", a);
        fieldCS.SetFloat("M", M);
        fieldCS.SetFloat("z", z);
        fieldCS.SetFloat("rim", rim);
        fieldCS.SetFloat("EHRadius", RAtZ(CalcEHRadius(M, a), z));
        fieldCS.SetBool("enableEventHorizon", enableEventHorizon);
        fieldCS.SetBool("enableErgosphere", enableErgosphere);
        fieldCS.Dispatch(
            kernel,
            Mathf.CeilToInt(width / 8f),
            Mathf.CeilToInt(height / 8f),
            1
        );

        GetComponent<Renderer>().material.mainTexture = target;
    }


    private void Update()
    {
        xpoint = xgrob + xfine;
        delT += Time.fixedDeltaTime;
        if ( Input.GetKeyDown(KeyCode.P) || (UpdateByTime && delT > UpdateTime)) //if delT = 1 for "continuous" updates 
        {
            texture = new Texture2D(width, height, TextureFormat.RGBA32, false);
            //texture.filterMode = FilterMode.Point;   // crisp pixels
            //texture.wrapMode = TextureWrapMode.Clamp;

            Dispatch();

            //GetComponent<Renderer>().material.mainTexture = texture;
            //delT = 0;
        }

    }

}

