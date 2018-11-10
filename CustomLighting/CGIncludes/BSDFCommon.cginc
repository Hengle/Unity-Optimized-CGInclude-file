#ifndef BSDFCOMMON
#define BSDFCOMMON
#include "MontcaloCommon.cginc"

//////////Diffuse BRDF
float Lambert(float NoL) {
	return NoL;
}

float LambertPi(float NoL) {
	return 0.3183098 * NoL;
}
#define POW5(x) x*x*x*x*x
float BurleyBRDF(float NoL, float NoV, float LoH, float roughness) {
    half FD90 = 0.5 + 2 * LoH * LoH * roughness;
    float NoL_Pow5 = POW5(1 - NoL);
    float NoV_Pow5 = POW5(1 - NoV);
    return (1 + (FD90 - 1) * NoL_Pow5) * (1 + (FD90 - 1) * NoV_Pow5) * NoL;
}
#undef POW5
//////////Cotton Cloth BRDF
float Inverse_GGX(float NoH, float Roughness) {
	float a2 = Roughness * Roughness;
	float A = 4;
	float d = (NoH - a2 * NoH) * NoH + a2;
	return rcp(1 + A * a2) * (1 + 4 * a2 * a2 / (d * d));
}

float Cloth_Geometry(float NoL, float NoV) {
	return rcp(4 * (NoL + NoV - NoL * NoV));
}

float3 Cloth_Cotton_BRDF(float NoH, float NoL, float NoV, float LoH, float roughness_Pow2, float3 specularColor) {
    float cloth_Cotton_GGX = Inverse_GGX(NoH, roughness_Pow2);
    float cloth_Cotton_Geometry = Cloth_Geometry(NoL, NoV);             
    float3 cloth_Cotton_Fersnel = FresnelTerm(specularColor, LoH);
    return (cloth_Cotton_GGX * cloth_Cotton_Geometry * cloth_Cotton_Fersnel) * NoL;
}

//////////Silk Cloth BRDF
void ConvertAnisotropyToRoughness(float roughness, float anisotropy, out float roughnessT, out float roughnessB) {
    float anisoAspect = sqrt(1 - 0.9 * anisotropy);
    roughnessT = roughness / anisoAspect; 
    roughnessB = roughness * anisoAspect; 
}

float3 ComputeGrainNormal(float3 grainDir, float3 V) {
	float3 B = cross(-V, grainDir);
	return cross(B, grainDir);
}

float3 GetAnisotropicModifiedNormal(float3 grainDir, float3 N, float3 V, float anisotropy) {
	float3 grainNormal = ComputeGrainNormal(grainDir, V);
	return normalize(lerp(N, grainNormal, anisotropy));
}

float TrowbridgeReitzAniso(float aniso, float roughness_Pow2, float NoH, float ToH, float BoH) {
    float aspect = sqrt(1 - aniso * 0.9);
    float X = max(0.001, roughness_Pow2 / aspect) * 5;
    float Y = max(0.001, roughness_Pow2 * aspect) * 5;
    return 1 / (3.1415926 * X*Y * Square(Square(ToH / X) + Square(BoH / Y) + NoH * NoH));
}

float Anisotropy_GGX(float ToH, float BoH, float NoH, float roughnessT, float roughnessB) {
    float f = ToH * ToH / (roughnessT * roughnessT) + BoH * BoH / (roughnessB * roughnessB) + NoH * NoH;
    return 1 / (roughnessT * roughnessB * f * f);
}

float Anisotropy_Geometry(float ToV, float BoV, float NoV, float ToL, float BoL, float NoL, float roughnessT, float roughnessB) {
    float aT = roughnessT;
    float aT2 = aT * aT;
    float aB = roughnessB;
    float aB2 = aB * aB;
    float lambdaV = NoL * sqrt(aT2 * ToV * ToV + aB2 * BoV * BoV + NoV * NoV);
    float lambdaL = NoV * sqrt(aT2 * ToL * ToL + aB2 * BoL * BoL + NoL * NoL);
    return 0.5 / (lambdaV + lambdaL);
}

float3 Cloth_Silk_BRDF(float NoL, float NoH, float LoH, float NoV, float ToH, float BoH, float ToV, float BoV, float ToL, float BoL, float roughnessT, float roughnessB, float3 specularColor) {
    float cloth_Silk_GGX = Anisotropy_GGX(ToH, BoH, NoH, roughnessT, roughnessB);
    float cloth_Silk_Geometry = Anisotropy_Geometry(ToV, BoV, NoV, ToL, BoL, NoL, roughnessT, roughnessB);             
    float3 cloth_Silk_Fersnel = FresnelTerm(specularColor, LoH);
    return (cloth_Silk_GGX * cloth_Silk_Geometry * cloth_Silk_Fersnel) * NoL;
}


//////////CookTorrance BRDF 
float Beckmann(float Roughness, float NoH) {
	float a = Roughness * Roughness;
	float a2 = a * a;
	float NoH2 = NoH * NoH;
	return exp((NoH2 - 1) / (a2 * NoH2)) / (PI * a2 * NoH2 * NoH2);
}

float GGX(float Roughness, float NoH) {
	float a = Roughness * Roughness;
	float a2 = a * a;
	float d = (NoH * a2 - NoH) * NoH + 1;	
	return a2 / (PI * d * d);					
}

float FresnelSchlickApprox(float HdV ,float F0){
    return F0 + (1 - F0) * pow(1 - HdV , 5);
}

float Vis_SmithJointApprox(float Roughness, float NoV, float NoL) {
	float a = Square(Roughness);
	float Vis_SmithV = NoL * (NoV * ( 1 - a ) + a);
	float Vis_SmithL = NoV * (NoL * ( 1 - a ) + a);
	return 0.5 * rcp(Vis_SmithV + Vis_SmithL);
}

float3 CookTorranceBRDF(float NoH, float NoL, float NoV, float LoH, float roughness_Pow2, float3 specularColor) {
    float pbr_GGX = GGXTerm(NoH, roughness_Pow2);
    float pbr_Geometry = SmithJointGGXVisibilityTerm(NoL, NoV, roughness_Pow2);             
    float3 pbr_Fersnel = FresnelTerm(specularColor, LoH);
    return (((pbr_Geometry * pbr_GGX) * UNITY_PI) * pbr_Fersnel) * NoL;
}


//////////Enviornment BRDF
float3 IntegrateBRDF(float Roughness, float NoV) {
	float3 V;
	V.x = sqrt( 1.0f - NoV * NoV );	
	V.y = 0;
	V.z = NoV;						

	float A = 0;
	float B = 0;
	float C = 0;

	const uint NumSamples = 64;
	[loop]
	for(uint i = 0; i < NumSamples; i++) {
		float2 E = Hammersley(i, NumSamples); 
		{
			float3 H = ImportanceSampleGGX(E, Roughness).xyz;
			float3 L = 2 * dot(V, H) * H - V;

			float NoL = saturate(L.z);
			float NoH = saturate(H.z);
			float VoH = saturate(dot(V, H));

			if(NoL > 0) {
				float Vis = Vis_SmithJointApprox(Roughness, NoV, NoL);

				float a = Square(Roughness);
				float a2 = a * a;
				float Vis_SmithV = NoL * sqrt( NoV * (NoV - NoV * a2) + a2 );
				float Vis_SmithL = NoV * sqrt( NoL * (NoL - NoL * a2) + a2 );
				float NoL_Vis_PDF = NoL * Vis * (4 * VoH / NoH);

				float Fc = pow(1 - VoH, 5);
				A += (1 - Fc) * NoL_Vis_PDF;
				B += Fc * NoL_Vis_PDF;
			}
		}

		{
			float3 L = CosineSampleHemisphere(E).xyz;
			float3 H = normalize(V + L);

			float NoL = saturate(L.z);
			float NoH = saturate(H.z);
			float VoH = saturate(dot( V, H ));

			float FD90 = (0.5 + 2 * VoH * VoH) * Roughness;
			float FdV = 1 + (FD90 - 1) * pow(1 - NoV, 5);
			float FdL = 1 + (FD90 - 1) * pow(1 - NoL, 5);
			C += FdV * FdL * (1 - 0.3333 * Roughness);
		}
	}
	return float3(A, B, C) / NumSamples;
}

float3 IntegrateBRDF_Cloth(float Roughness, float NoV) {
	float3 V;
	V.x = sqrt( 1.0f - NoV * NoV );	
	V.y = 0;
	V.z = NoV;						

	float A = 0;
	float B = 0;
	float C = 0;

	const uint NumSamples = 64;
	[loop]
	for(uint i = 0; i < NumSamples; i++) {
		float2 E = Hammersley(i, NumSamples); 
		{
			float3 H = ImportanceSampleInverseGGX(E, Roughness).xyz;
			float3 L = 2 * dot(V, H) * H - V;

			float NoL = saturate(L.z);
			float NoH = saturate(H.z);
			float VoH = saturate(dot(V, H));

			if(NoL > 0) {
				float Vis = Cloth_Geometry(NoV, NoL);

				float a = Square(Roughness);
				float a2 = a * a;
				float Vis_SmithV = NoL * sqrt(NoV * (NoV - NoV * a2) + a2);
				float Vis_SmithL = NoV * sqrt(NoL * (NoL - NoL * a2) + a2);
				float NoL_Vis_PDF = NoL * Vis * (4 * VoH / NoH);

				float Fc = pow(1 - VoH, 5);
				A += (1 - Fc) * NoL_Vis_PDF;
				B += Fc * NoL_Vis_PDF;
			}
		}

		{
			float3 L = CosineSampleHemisphere(E).xyz;
			float3 H = normalize(V + L);

			float NoL = saturate(L.z);
			float NoH = saturate(H.z);
			float VoH = saturate(dot(V, H));

			float FD90 = ( 0.5 + 2 * VoH * VoH ) * Roughness;
			float FdV = 1 + (FD90 - 1) * pow(1 - NoV, 5);
			float FdL = 1 + (FD90 - 1) * pow(1 - NoL, 5);
			C += FdV * FdL * (1 - 0.3333 * Roughness);
		}
	}
	return float3(A, B, C) / NumSamples;
}

float3 PreintegratedGF_PC(float3 SpecularColor, float Roughness, float NoV) {
	float2 AB = IntegrateBRDF_Cloth(Roughness, NoV).xy;
	return SpecularColor * AB.x + saturate(50 * SpecularColor.g) * AB.y;
}

float3 PreintegratedGF_PC_Cloth(float3 SpecularColor, float Roughness, float NoV) {
	float2 AB = IntegrateBRDF(Roughness, NoV).xy;
	return SpecularColor * AB.x + saturate(50 * SpecularColor.g) * AB.y;
}

float3 PreintegratedGF_LUT(sampler2D PreintegratedLUT, float3 SpecularColor, float Roughness, float NoV) {
	float2 AB = tex2Dlod(PreintegratedLUT, float4(clamp(Roughness, 0.001, 0.999), NoV, 0, 0));
	return SpecularColor * AB.x + saturate(50 * SpecularColor.g) * AB.y;
}

half3 PreintegratedGF_Mobile(half3 SpecularColor, half Roughness, half NoV) {
    const half4 c0 = {-1, -0.0275, -0.572, 0.022};
    const half4 c1 = {1, 0.0425, 1.04, -0.04};
    half4 r = Roughness * c0 + c1;
    half a004 = min(r.x * r.x, exp2(-9.28 * NoV)) * r.x + r.y;
    half2 AB = half2(-1.04, 1.04) * a004 + r.zw;
    AB.y *= saturate(50 * SpecularColor.g);
    return SpecularColor * AB.x + AB.y;
}

half PreintegratedGF_Mobile_Nonmetal(half Roughness, half NoV) {
	const half2 c0 = { -1, -0.0275 };
	const half2 c1 = { 1, 0.0425 };
	half2 r = Roughness * c0 + c1;
	return min( r.x * r.x, exp2( -9.28 * NoV ) ) * r.x + r.y;
}
#endif