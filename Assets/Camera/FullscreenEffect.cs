using UnityEngine;

[ExecuteInEditMode]
public class FullscreenEffect : MonoBehaviour
{
    public Material material;

    void OnRenderImage(RenderTexture src, RenderTexture dest)
    {
        if (material != null)
            Graphics.Blit(src, dest, material);
        else
            Graphics.Blit(src, dest);
    }

    void Update()
    {
        material.SetFloat("_Zoom", Mathf.Pow(1.01f, Time.time * 10)); // zoom in over time
    }

}
