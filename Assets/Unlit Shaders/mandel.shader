Shader "Custom/mandel"
{
    Properties
    {
        _Offset("Offset", Vector) = (0,0,0,0)
        _Zoom("Zoom", Float) = 500
        _MaxIterations("Max Iterations", Float) = 1000

    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Geometry" }
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

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            float _Zoom;
            float4 _Offset;
            float _MaxIterations;

            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float2 c;
                c.x = (i.uv.x - 0.5) * 3.5 / _Zoom + _Offset.x;
                c.y = (i.uv.y - 0.5) * 2.0 / _Zoom + _Offset.y;

                float2 z = float2(0,0);
                int iter = 0;

                for(int n = 0; n < _MaxIterations; n++)
                {
                    float x = z.x*z.x - z.y*z.y + c.x;
                    float y = 2*z.x*z.y + c.y;
                    z = float2(x,y);
                    if(dot(z,z) > 4) break;
                    iter++;
                }

                float t = iter / _MaxIterations;

                // Smooth coloring: map t to blue-white gradient
                float3 color = lerp(float3(0,0,0.2), float3(1,1,1), sqrt(t));
                return fixed4(color,1);
            }

            ENDCG
        }
    }
}
