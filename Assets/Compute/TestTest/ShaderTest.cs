using UnityEngine;

public class ShaderTest : MonoBehaviour
{
    public ComputeShader CS;
    RenderTexture target;

    void Update()
    {
        if (Input.GetKeyDown(KeyCode.P))
        {
            target = new RenderTexture(256, 256, 0, RenderTextureFormat.ARGBFloat);
            target.enableRandomWrite = true;
            target.Create();

            int kernel = CS.FindKernel("CSMain");
            Debug.Log("Kernel index: " + kernel);

            CS.SetTexture(kernel, "Result", target);
            CS.Dispatch(kernel, 256 / 8, 256 / 8, 1);

            Debug.Log("Done. target created: " + target.IsCreated());
        }
    }

    void OnGUI()
    {
        if (target != null)
            GUI.DrawTexture(new Rect(0, 0, 256, 256), target);
    }
}
