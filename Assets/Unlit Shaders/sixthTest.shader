Shader "Custom/FullscreenTest"
{
    Properties
    {
        _Color("Color", Color) = (1,0,0,1)
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct circle
            {
               float2 center;
               float radius;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            fixed4 _Color;

            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                if((i.uv.x-0.5) * i.uv.x + i.uv.y * i.uv.y < 0.25)
                {
                    return fixed4(1, 0, 0, 1); // Red color inside the circle
                }
                return fixed4(0, 0, 0, 1);
            }
            ENDCG
        }
    }
}