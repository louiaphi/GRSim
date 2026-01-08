using UnityEngine;

public class NewMonoBehaviourScript : MonoBehaviour
{
    public ComputeShader cs;
    RenderTexture rt;

    void Start()
    {
        rt = new RenderTexture(512, 512, 0);
        rt.enableRandomWrite = true;
        rt.Create();

        int kernel = cs.FindKernel("CSMain");
        cs.SetTexture(kernel, "Result", rt);
        cs.Dispatch(kernel, 512 / 8, 512 / 8, 1);

        GetComponent<Renderer>().material.mainTexture = rt;
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
