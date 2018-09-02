# Optimized CGInclude files
## Introduction:
### Standard Lighting Model:
* In Unity, The Standard PBR Lighting Model use GGX algorithm as direct light specular and surfaceproduction simulate as indirect reflection specular. However, it looks "tough" and not "comfortable" for some artists who used to talked about that with me. Thus, I determine to create this project to transform Unity's standard lighting model into custom lighting model.
* For direct light, an optimized Cook Torrence algorithm has been provided, which take less instructions than the official one with the same (or similar) result. For indirect light, we use Pre-Integrate ramp texture to simulate the surface production and attenuation by using Montcalo Intergration.
![demo]("demo1.png")
![demo]("demo2.png")
### PCSS Directional Shadow:
* We use NVIDIA PCSS to replace Unity's soft shadow algorithm. By randomly sample 16, 32 or 64 times, the edge of shadow will looks much better. We strongly suggest that users should add Temporal AA post processing effects in the scene. The Temporal Filter will reduce the noise of the soft shadow.
![demo]("demo0.png")
## Manual:
### Standard Lighting Model:
1. Copy and replace the files in the folder "CGIncludes" to Unity/Editor/Data/CGIncludes(Remember to make backup);
2. Drag "Resources" and "PreIntSpecular.cs" into your scene and enable the component in a random GameObject, usually the "PreIntSpecular.cs" should be singleton in the whole game.
3. The Lighting Model should support Standard shaders and regular surface shaders. If you are writing your own VF PBR shaders, make sure you are calling the methods in UnityStandardBRDF.cginc.
4. Enjoy your new lighting!

### PCSS Directional Shadow
1. Drag the folder "Directional Shadow" into your Unity project. 
2. Open "Edit/Project Settings/Graphics" interface.
3. Replace the "Screen space shadows" to custom shader and add "Hidden/PCSS_Directional" as the custom shader.
4. Enable the component "PCSS_Directional.cs" in your scene. This component should be singleton in one scene(not in the whole game).

## Supported Version:
* This plugin should support Unity 2017, Unity 2018 and newer version.