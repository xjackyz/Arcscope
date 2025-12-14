ArcScope 1.0 设计文档

目录
•	0. 项目概述
•	1. 自适应总原则
•	2. 系统架构与核心功能
•	3. 分析与归一化流程
•	4. 曲线体系
•	5. 多尺度结构分段
•	6. 自动诊断
•	7. FilmEngine 设计
•	8. 数据与持久化（SQLite）
•	9. Bridge 层
•	10. macOS App
•	11. 开发优先级

0. 项目概述

项目名称：ArcScope
目标平台：macOS（Apple Silicon 优先）

核心定位：
•	面向电影创作者（导演、剪辑、制片、研究者）的本地电影结构分析工具。
•	将一部影片转化为若干条随时间变化的"结构曲线"，帮助创作者看清节奏、情绪、信息和视听关系。
•	尤其用于：
    o	剪辑/混音阶段的复盘与调节
    o	对经典影片的结构"拆解学习"
    o	自己作品的二刷分析与问题段定位

关键特征：
•	完全本地分析（无云端上传，规避版权问题）
•	统一时间轴：所有曲线以秒为横坐标、[0,1] 归一化为纵坐标
•	多模态综合：剪辑节奏、运动、声音、色彩、字幕信息、表演情绪
•	自动诊断：标出"可能拖沓 / 情绪空转 / 视听错位 / 解释压死情绪"等区段
•	参考对比：可叠加经典影片曲线做结构对照

产品策略：
•	ArcScope 1.0 只发布 macOS 版本
•	UI 完全使用 SwiftUI，提供原生 macOS 交互体验
•	核心分析引擎采用平台无关的 C++ 设计，便于后续扩展

架构策略：
•	统一视频处理流：不区分 ProRes / H.264 / HEVC 等编码格式，所有影片解码后统一下采样到固定高度（默认 720p），使用固定分析帧率（默认 10 fps）进行所有特征提取，最终在 1 Hz 上输出曲线。保证结果可比性、性能可控，避免高码率格式导致的分析爆炸
•	稳健自适应归一化：所有特征均使用中位数 + MAD + Sigmoid 的统一公式映射到 [0,1]，无需手动阈值，自动适配每部影片
•	平台解耦设计：核心分析引擎（FilmEngine）与平台特定代码完全解耦，便于未来移植到命令行工具或其他平台

注：ArcScope 明确区分 Core Curves 与 Metadata Overlays：前者（Pace/Sound/Color/Info/Arousal/Face）参与 PCA 与 Diagnostics；后者（BPM、Key、ShotScale、CameraMotion、IsSpeaking 等）仅作可视化注释，从不进入主算法。

1. 自适应总原则
1.1 ArcScope 自适应金律
•	禁止写死绝对阈值、固定权重、固定区间长度，所有“高低 / 强弱 / 突变”判断必须从影片自身分布推导
•	统一遵循“statistically adaptive + 互相归一 + 互相参照 + 全流程无硬编码”的策略，不允许 θ=0.7、T_min=30 秒、权重 0.2/0.5 之类的魔法数
•	曲线、Issue、平滑窗口、诊断持续时间等都交给 robust statistics + PCA 等统计方法来确定，影片风格自己定义标尺

1.2 实施要点
•	稳健中位数 + MAD + Sigmoid 提供所有特征的无参数归一化
•	多模态特征统一做 PCA 投影，取主方向而非手调权重
•	诊断阈值与持续时间基于分位数 + 稳定段长度，任何“30 秒”“0.7”式常数都禁止出现
•	窗口长度、平滑系数等“时间”参数均由影片节奏 λ 自动映射，确保快片/慢片都能匹配
•	情绪分析遵循叙事主视角：只跟踪“主角”这条 dominant face trajectory，其他人脸全部视为噪声

2. 系统架构与核心功能
2.1 分层架构
ArcScope 由五层组成：
1.	macOS App 层（SwiftUI）
o	影片导入、分析控制
o	曲线可视化与问题段标记
o	参考影片选择与对比
o	报告视图
2.	Bridge 层（Objective-C++）
o	Swift ↔ C/C++ 调用桥接
o	封装 Apple 专有框架：
	AVFoundation（解码/音频）
	Vision（光流、检测人脸）
	CoreML（情绪/表情模型）
	Accelerate / vDSP（FFT、滤波）
3.	分析引擎层（C++17，FilmEngine）
o	媒体读取与特征抽取（与 Bridge 协作）
o	各类低层特征：剪辑、运动、声音、色彩、字幕信息
o	曲线构造（Pace / Sound / Color / Info / Arousal / FaceAffect）
o	自动诊断引擎（IssueSegment）
4.	本地数据库层（SQLite）
o	films：影片信息
o	curves：曲线数据
o	metrics：标量指标（均值、方差、峰值统计等）
o	segments：问题区段诊断结果
o	references_library：用户勾选的“参考影片”
5.	模型与参考库层
o	CoreML 模型（表情/情绪、可选多模态 arousal 模型）
o	用户自行生成的“经典影片曲线数据库”
2.2 数据流
1.	用户在 SwiftUI 中选择本地影片文件。
2.	Swift 调用 Bridge C 接口 arcscope_analyze_film。
3.	Bridge 使用 AVFoundation + Vision + CoreML 提取：
o	每秒光流强度、帧图像缓冲、音频 PCM 片段、人脸情绪特征等
4.	C++ FilmEngine 通过 MediaReader/特征模块构造各种 Feature(t)：
o	Shot 时长、cut 密度、运动强度
o	音频 RMS / 频带能量
o	色彩亮度/饱和度/暖度、灰度感
o	字幕文本字数 / exposition 比例
o	Face presence / 表情 intensity / valence / arousal
5.	CurveEngine 使用统一公式将这些特征组合成 5–6 条归一化曲线。
6.	DiagnosticsEngine 在曲线空间中扫描，找出若干 IssueSegment。
7.	FilmEngine 将结果写入 SQLite（films/curves/metrics/segments）。
8.	SwiftUI 读取数据库，显示曲线、问题段，并支持与参考影片曲线叠加。
 
2.3 功能清单
2.3.1 核心功能（1.0 必须实现）
•	影片导入与分析启动
•	5 条核心曲线：
1.	Pace(t)：节奏 / 剪辑 + 运动
2.	SoundEnergy(t)：声音能量 / 压力
3.	ColorMood(t)：色彩情绪 / 视觉温度
4.	InfoDensity(t)：信息密度 / 说明性负荷
5.	Arousal(t)：综合唤醒曲线（多模态）
•	第 6 条核心曲线：FaceAffect(t)（表演存在度与情绪强度）
◦	当检测到稳定的主角轨迹时启用
◦	若整片无可靠主角轨迹，则自动禁用：
‣	数据库不写入 face_affect 曲线
‣	Arousal 不使用面部特征
‣	UI 中隐藏 Face 曲线开关，并在信息栏提示 “No stable dominant face; arousal driven by pace+sound+info”
•	自动诊断 4 类区段：
o	低活跃盆地（LowActivity）
o	节奏过猛但情绪空转（OvercutFlat）
o	视听错位（AudioVisualMisalign）
o	解释压死情绪（InfoKillsEmotion）
•	曲线总览 UI：
o	统一时间轴，0–1 纵轴，曲线可勾选显示
o	支持缩放 / 拖动
•	问题区段列表 + 时间轴高亮
•	本地影片库（已分析影片列表）
•	将任意影片标记为 “参考影片”，在曲线视图中叠加显示
2.3.2 高级功能（1.0 底座预留接口，可后续增强）
•	使用 CoreML 多模态 arousal 模型替代手工加权部分
•	Dialog 情绪细分（文本 valence/arousal）
•	类型/导演风格地图（PCA/UMAP 空间可视化）
•	报告导出（PDF / Markdown）
 
3. 分析与归一化流程
3.1 统一视频处理策略
ArcScope 对所有视频格式采用统一的分析流程，从产品角度只提供一个"完整模式"，避免用户面对多个版本选择。

3.1.1 视频解码与下采样
•	不区分编码格式：ProRes、H.264、HEVC 等统一处理
•	解码后立即下采样：
o	固定高度：720p 或 1080p（可配置，默认 720p）
o	保持原始宽高比
•	固定分析帧率 f_analysis：
o	默认 10 fps（可配置范围 5–30 fps）
o	所有内部特征分析均在此帧率上进行

3.1.2 各类特征的采样策略
基于 f_analysis = 10 fps 的内部流：

1.	光流分析
o	在 10 fps 连续帧上计算光流
o	每秒得到 10 个 flow_mag 值
o	汇总为该秒的平均运动强度 M(k)

2.	表情分析（降频 + gating）
o	人脸/表情推理采用降频策略：
    	基础 f_analysis = 10 fps（光流/色彩等）
    	f_face = \min(3\,\text{fps}, f_{analysis}/stride_{face})，stride_{face} ∈ {3,4}，即每 3–4 帧才送一次检测 + 表情
    	每秒最多挑选 3 帧进入 CoreML 表情模型
o	结合运动 gating：若该秒的光流均值 M(k) < θ_M（θ_M = Q_M(0.6)）且主角轨迹稳定，则复用上一秒的表情结果；只有当 M(k) ≥ θ_M 或 FacePresence 出现 0↔1 切换时才重新推理
o	汇总为该秒的 face presence、expression intensity、valence/arousal

3.	色彩分析
o	每秒从 10 帧中均匀抽取 1–2 帧
o	统计亮度、饱和度、暖度等特征

4.	音频分析
o	按原始音频采样率处理
o	以 1 秒为单位计算 RMS / 频带能量

5.	字幕/剪辑分析
o	基于原始时间戳
o	按 1 秒窗口汇总字数、镜头边界等

3.1.3 输出时间轴
•	所有曲线最终按 1 Hz 采样输出：
[
t_k = k,\quad k = 0,1,\dots, N-1
]
•	默认将第 k 秒视为区间 [k,k+1]，其中心时间为 t_k^{center} = k + 0.5；所有 1 Hz 特征（Pace/Sound/Face 等）均代表该中心时间
•	每一秒的 motion / face / sound / color 等数据均从内部 10 fps 分析结果汇总而来，并在中心时间对齐
•	用户看到的始终是统一的 1 秒采样曲线，无需关心内部分析帧率

工程优势：
•	自动降分辨率 + 降帧率，无需针对不同编码格式设置不同分支
•	分析性能可控，避免 ProRes 等高码率格式导致的性能爆炸
•	所有影片曲线具有可比性

3.2 归一化方式：中位数+MAD+Sigmoid（唯一方案）
ArcScope 所有连续特征统一采用以下公式归一化到 (0,1)：
对每个序列 X_k：
[
m = \text{median}(X_k), \quad s = \text{median}(|X_k - m|)
]
[
Z(k) = \frac{X_k - m}{s + 10^{-6}}, \quad \tilde{X}(k) = \frac{1}{1 + e^{-Z(k)}}
]
此方法兼具稳健性（中位数/MAD）与平滑压缩（Sigmoid），无分支、无可选项。

归一化的含义：
•	不表达"绝对高低"，而是"在这部片自己的分布里，这一秒是高是低"
•	整片平淡的影片：曲线接近一条 0.5 的平线
•	局部有相对明显变化：Z 提起，但不会被吹到天上
•	用户无需调节参数，算法自动适配每部片的风格

最终结果：
所有曲线呈现为 [0,1] 范围，便于叠加对比，且对不同风格影片自适应。

3.3 节奏与情绪动力学基础
3.3.1 自适应平滑窗口（统一公式）
所有基于时间窗口的局部统计（Pace 局部 ASL / CutDensity、多尺度变点检测等）使用同一公式确定窗口长度：
1.	取 Pace 的稳健 z 序列 Z_Pace(t)
2.	计算平均斜率：
[
\lambda = \mathbb{E}_t \big[ |Z_{Pace}(t) - Z_{Pace}(t-1)| \big]
]
3.	映射并截断：
[
W^* = 12 e^{-2\lambda}, \quad W = \min(\max(W^*, 4), 30) \text{ 秒}
]
快片（λ 大）→ W 变短；慢片 → W 变长。文档中出现的 W/Δ 均指此窗口。

3.3.2 情绪动力学指标
ArcScope 后续的结构分析会使用 Arousal(t) 的变化率：
•	一阶导（情绪斜率）
[\
A'(k) = A(k) - A(k-1)
]
•	二阶导（曲率 / 加速度）
[\
A''(k) = A(k+1) - 2A(k) + A(k-1)
]
•	节奏-情绪分歧：
[\
D_{PA}(k) = Z_{Pace}(k) - Z_{Arousal}(k)
]
这些量在 Scene/Sequence segmentation、诊断等模块中用于检测情绪反转、build-up/释放、节奏失配等结构点。

3.4 色彩处理基线（Color Science Pipeline）
ArcScope 的色彩特征遵循业内标准的“双层”流程：
1.	色度学统一：所有视频帧先转换到场景线性的 ACEScg（或同级 Rec.2020 线性）空间，确保不同相机/色域在同一基准下比较。
2.	感知空间投影：再将线性 RGB 映射到感知均匀模型（默认 CAM16-UCS，可选 OKLab/ICtCp）。计算以下量：
    •	感知亮度 J（或 L*）
    •	色度/饱和度 M、C
    •	色调角 h
    •	色差 ΔE（镜头间颜色跳变）
所有颜色统计（均值、方差、色调分布、色差轨迹）都在该感知空间里计算，保证“数值变化” = “人眼感知变化”。色彩和情绪、和谐/紧张度等指标也基于 CAM16/OKLab 的度量来定义。

3.4.1 理论先验 + 片内校准（Color-in-Context）
•	ArcScope 提供一组默认的 color→emotion 映射：基于亮度/饱和度/暖度推断 Valence_color^{raw}(t)、Arousal_color^{raw}(t)
•	为适应具体影片的语境，系统允许 per-film calibration：用户可在时间轴上标注少量情绪片段（valence/arousal），算法对 ColorWarmth/ColorEnergy 等特征做线性回归/小模型校正，得到 Valence_color(t)、Arousal_color(t) 的片内定制版本
•	所有 color→emotion 关系都遵循 Color-in-Context Theory：默认先验 + 片内 calibration = 最终情绪特征

3.4.2 Psychoacoustic Loudness & Clarity（声音预处理）
声音部分的输入也引入人耳模型：
•	主观响度 Loudness(t)：对频谱应用 A-weighting/RLB weighting，按 ITU-R 1770/LUFS 思路整合，形成更贴近人耳的声压特征
•	对白清晰度 DialogueClarity(t)：结合语音频段能量（300–3400 Hz）与背景噪声能量，计算 SNR；可与 SpeechDuty 共同用于 Info/Sound
•	SpectralBrightness(t)：使用谱重心/高频能量比，反映声音的“亮/刺耳”程度
这些特征在 4.2 中作为子曲线的一部分，保证 SoundEnergy 既有物理也有人耳感知依据。

3.4.3 Subtitle & NLP Pipeline（信息预处理）
ArcScope 的 InfoDensity 依赖离线文本处理：
1.	字幕/ASR 解析：将字幕或语音识别结果解析为 SubtitleSegment(start_sec, end_sec, text)
2.	NLP 预处理：使用 macOS `NaturalLanguage`（或可替换的 CLI）对每个段落执行：
    •	语言识别 → 若为英文/欧语，使用 NLTagger 提供的词性/lemma；若为中文等，至少提供分词
    •	分词、词性、lemma、NER 标记，构建 Token(form, lemma, POS, flags)
    •	句子划分，生成 Sentence(t_start, t_end, tokens)
3.	每秒聚合：将 Token/Sentence 投射到秒级，计算 Words(k)、ConceptDensity(k)、NewConceptRatio(k)、SentenceComplexity(k)、ExpoScore(k) 等特征
4.	所有特征经中位数+MAD 归一化后进入 VerbalLoad/InfoDensity 的 PCA
该 pipeline 保持后端可替换（Apple NL、Python CLI 等），输出接口固定，确保 ArcScope 可以扩展多语种和更高级别的语言特征。
 
4. 曲线系统设计
下面每条曲线都包含：
•	研究和直觉来源
•	影视语境下的意义
•	精确算法定义（可直接实现）
4.1 节奏 / 剪辑曲线 Pace(t)
4.1.1 含义
综合反映 剪辑速度 + 镜头运动 的“节奏感”：
•	镜头越短、切得越密、画面运动越强 → Pace 越高。
4.1.2 输入特征
•	镜头集合：
[
\mathcal{S} = {(t_{s,i}, t_{e,i})}_{i=1}^M
]
•	每个镜头时长：
[
d_i = t_{e,i} - t_{s,i}
]
•	每帧光流强度 (flow_mag(f))（由 Vision 光流计算）
窗口宽度 (\Delta) 统一使用 3.3 节定义的 W（由 Pace z 序列斜率得出，截断在 4–30 秒）。
4.1.3 局部特征
1.	局部 ASL（Average Shot Length）
对于时间 (t)，窗口
[
W(t) = [t - \Delta/2,, t + \Delta/2]
]
内的镜头索引集合：
[
S(t) = { i \mid [t_{s,i}, t_{e,i}] \cap W(t) \neq \varnothing }
]
则：
[
ASL(t) = \frac{1}{|S(t)|} \sum_{i\in S(t)} d_i
]
2.	剪辑密度 CutDensity
[
CD(t) = \frac{#{\text{cut in } W(t)}}{\Delta}
]
3.	运动强度 Motion
令 (F(t)) 为窗口内的帧集合：
[
M(t) = \frac{1}{|F(t)|} \sum_{f\in F(t)} flow_mag(f)
]
对上述三个特征分别做稳健 z-score（按第 3.2 节方法：中位数 + MAD），得到 (Z_{ASL}, Z_{CD}, Z_{M})。
4.1.4 曲线定义：主方向自适应
ArcScope 不再手写 (\alpha) 权重，而是在 [ -Z_{ASL}, Z_{CD}, Z_{M} ] 的多维轨迹上做 PCA，自动找出“节奏主方向”。

1.	构造矩阵 (\mathbf{F}(t) = [-Z_{ASL}(t), Z_{CD}(t), Z_{M}(t)])
2.	对 \mathbf{F} 的列再次稳健标准化，计算协方差矩阵
3.	取第一主成分 (w_{pace})
4.	若 (w_{pace}) 导致 Pace 趋势与直觉反向，则整体乘以 (-1)，保持“越快越高”

[
\text{Pace_raw}(t) = \mathbf{F}(t) \cdot w_{pace}
]
再对 (\text{Pace_raw}) 做一次稳健 z-score + sigmoid（按第 3.2 节），得到：
[
Pace(t) \in (0,1)
]
 
4.2 声音曲线体系
4.2.1 Loudness(t)
•	依据 3.4.2 的 psychoacoustic loudness（A/RLB weighting + LUFS 风格整合），得到 Z_{Loudness}(k)，归一化成 Loudness(t)，代表人耳感知的声压/压迫感。

4.2.2 RhythmDrive(t)
•	利用谱通量 SF[m] 与主节拍附近能量占比计算 RhythmEnergy(k)，经归一化后得到 RhythmDrive(t)，刻画节奏驱动力。

4.2.3 HarmonicTension(t)
•	基于 Chroma 与当前 key triad 的偏离：
[\
\text{Tension}(k) = 1 - \frac{\sum_{p\in P_{triad}(K^*(k))} c_k[p]}{\sum_p c_k[p] + \varepsilon}
]
归一化后得到 HarmonicTension(t)，反映音乐紧张度。

4.2.4 SpectralBrightness(t)
•	使用谱重心或高频能量比（Z_{HF}) 表示声音的明亮/刺耳程度，归一化得到 SpectralBrightness(t)。

4.2.5 综合曲线：SoundEnergy(t)
构造向量：
[\
\mathbf{S}(k) = [Z_{Loudness}(k), Z_{RhythmDrive}(k), Z_{Tension}(k), Z_{Brightness}(k)]
]
做 PCA，取第一主成分 (w_{sound})：
[\
\text{Sound_raw}(k) = \mathbf{S}(k) \cdot w_{sound}
]
若与 Z_{Loudness} 相关为负则乘 (-1)。归一化得到 SoundEnergy(t_k)，作为 Arousal/Diagnostics 的声音主轴。

4.2.6 声部角色（metadata / 未来特征）
•	DialogueEnergy(k)、MusicEnergy(k)、EffectsEnergy(k)：按声部分离（对白频段能量、音乐通道、残余 = Effects）
•	DialogueDominance(t) = DialogueEnergy / (总能量 + ε)，Music/Effects 同理
•	这些曲线可作为 metadata overlay，或在诊断中用于解释 InfoKillsEmotion/AV Misalign 的来源

Metadata（仅显示，不入主轴）：BPM(t)、TempoConfidence(t)、Key(t)/KeyChangeDensity(t)、Dialogue/Music/Effects Dominance 条带。
 
4.3 色彩氛围与和谐度曲线

4.3.1 ColorWarmth(t)
描述画面暖冷偏移：在 CAM16/OKLab 中统计暖色像素占比、色调角分布，计算 Z_{Warm}(k)，再按 3.2 归一化得到 ColorWarmth(t)。

4.3.2 ColorEnergy(t)
由感知亮度与饱和度联合反映“视觉能量”：构造 [Z_{J}(k), Z_{Sat}(k)]，取 PCA PC1 作为 ColorEnergy_raw(k)，归一化成 ColorEnergy(t)。用于表示“亮且鲜 vs 暗且灰”的节奏。

4.3.3 ColorHarmonyScore(t)（metadata）
在 CAM16-UCS 中对每秒的像素做 K-medians 聚类，计算色群均值间的距离与群内方差，合成 [0,1] 的和谐度：群间差距越大、群内越集中 → 和谐度高；色调冲突越大 → 和谐度低。默认作为 metadata overlay，也可在未来升级为曲线。

4.3.4 Color–Affect Binding
依据 color psychology，ArcScope 为每部片自适应地学习颜色对 valence/arousal 的贡献：
•	暖度/亮度 → valence 预测：
[\
\hat{V}_{color}(k) = a_w W(k) + a_b B(k)
]
•	饱和度 → arousal 预测：
[\
\hat{A}_{color}(k) = a_c C(k)
]
系数 a_w,a_b,a_c 通过对 Scene 平均值做线性回归或 PCA 获得，完全自适应。颜色与真实情绪的对齐度：
[\
Align(k) = \big(Z_{FaceValence}(k) - Z_{\hat{V}_{color}}(k)\big)^2 + \big(Z_{Arousal}(k) - Z_{\hat{A}_{color}}(k)\big)^2
]
归一化后得到 AlignNorm(k)，越小表示颜色与情绪一致，越大表示 counterpoint 紧张。ColorMood(t) 最终可视为 Z_{Warm},Z_{J},Z_{Sat},Z_{Gray},HarmonyScore,Align 等特征的 PC1，供 Arousal/Diagnostics 使用。

4.3.5 ColorState(t) 与模式切换
在 CAM16/OKLab 的 [Warmth,Brightness,Saturation] 空间对全片帧做 KMeans (K=3–5)，得到离散色彩模式 State(t)。
•	StateShift(t) = 1[State(t) ≠ State(t-1)] 用于检测 palette 切换
•	可在 Scene/Sequence 指纹中存储每段的主色模式、色彩指纹（Hue 直方图、ColorWarmth/Energy/Harmony 均值）

Metadata（可选）：Hue track、ColorHarmonyScore、Character Palette、ColorState 色带、AlignNorm(k) 曲线、Counterpoint 区域高亮等。
 
4.4 信息密度曲线 InfoDensity(t)
4.4.1 含义
衡量每秒输入的语言/视觉/事件量，反映认知负荷；信息过载会挤压情绪/节奏，尤其在 Arousal 低时。

4.4.2 VerbalLoad(t)
语言负荷向量：
[\
\mathbf{X}_{verbal}(k) = [Z_{Words}, Z_{SpeechDuty}, Z_{ConceptDensity}, Z_{NewConceptRatio}, Z_{SentenceComplexity}, Z_{Exposition}, Z_{DialogueClarity}]
]
这些特征来自 3.4.3 的 NLP pipeline。PCA → VerbalLoad_raw(k)，归一化 → VerbalLoad(t)。

4.4.3 VisualLoad(t)
[\
\mathbf{X}_{visual}(k) = [Z_{VisualEntropy}, Z_{ShotScaleNumeric}, Z_{CameraMotionComplexity}, Z_{FaceChange}]
]
PCA → VisualLoad(t)，可解释为“画面复杂度/读取成本”。

4.4.4 EventLoad(t)
[\
\mathbf{X}_{event}(k) = [Z_{SoundEvents}, Z_{MusicChange}, Z_{CutDensity}, Z_{ObjectMotionJump}]
]
PCA → EventLoad(t)。

4.4.5 ExpositionCurve(t)
ExpoScore 基于数字/专名/术语/句长/复杂度等特征打分（见 3.4.3 输出），归一化后形成 ExpositionCurve(t)，既可独立展示，也纳入 VerbalLoad/InfoDensity。

4.4.6 InfoDensity(t)
组合 Verbal/Visual/Event/Exposition：
[\
\mathbf{X}_{info}(k) = [Z_{VerbalLoad}, Z_{VisualLoad}, Z_{EventLoad}, Z_{Exposition}]
]
PCA PC1 → Info_raw(k)，若与 Z_{VerbalLoad} 相关 < 0 则乘 (-1)；归一化后得到 InfoDensity(t)。

4.4.7 Cognitive Difficulty（预留）
ArcScope 预留更高级的认知特征接口（ConceptDensity、NewConceptRatio、SentenceComplexity 已写入 4.4.2；未来可继续扩展 discourse/指代复杂度）。
 
4.5 表演存在度 / 情绪曲线 FaceAffect(t)
4.5.1 叙事视角：只看主角
ArcScope 分析的是叙事主视角，而非群体心理。电影语言（镜头距离、出镜时长、画面占比）早已告诉我们谁是主角，系统只跟踪这条 dominant face trajectory，其余人脸全部视为噪声：
•	吵架戏/群戏 → 依然只跟主角表情，整体张力交由 Pace/Sound/Motion 反映
•	主角缺席 → FaceAffect 自动降为 NULL/低权重，不让背景路人拉高情绪

4.5.2 主角轨迹识别（全自适应）
每帧检测所有人脸并进行多目标跟踪，得到若干轨迹 track_i。对每条轨迹统计：
•	出现时长比例 T_i = frames_i / total_frames
•	平均脸面积占比 S_i（脸像素 ÷ 画面像素）
•	检测置信度均值 C_i（来自人脸/表情模型）
定义信号质量：
[
\text{Score}_i = T_i \cdot S_i \cdot C_i
]
为了自适应地筛掉噪声轨迹：
1.	计算所有轨迹的 T_i 与 C_i 的 75% 分位点：T_q = Q_T(0.75)、C_q = Q_C(0.75)
2.	只有满足 T_i ≥ T_q 且 C_i ≥ C_q 的轨迹才进入候选集
3.	在候选集中取 Score_i 最大者作为主角轨迹 track_*
4.	若没有轨迹同时通过 T/C 分位门槛，则认定“无可靠主角轨迹”，整片关闭 FaceAffect：
    •	curves 表不写入 curve_type='face_affect'
    •	仅写入 curve_type='face_presence'，值全为 0
    •	Arousal 不包含 face 维度
    •	UI/报告文案：No stable dominant face detected；Arousal based on pace+sound+info

4.5.3 曲线与 mask：主角情绪强度
当主角轨迹在当前秒可见（且满足 3.1.2 中的 motion gating 条件）时：
•	抽取主角的 face presence P_*(k)、expression intensity E_*(k)、模型 arousal A_*(k)（valence 另行输出）、脸面积占比 S_*(k)
•	分别计算 Z_P(k)、Z_E(k)、Z_A(k)（按 3.2 归一化）
•	定义景别权重：
[
w_{close}(k) = \text{clip}\Big( \frac{S_*(k)}{\text{median}(S_*) + 10^{-6}}, 0, 1 \Big)
]
•	原始面部信号：
[
\text{Face\_raw}(k) = w_{close}(k) \cdot (Z_P(k) + Z_E(k) + Z_A(k))
]
•	FacePresence(k) = 1

当主角不在场或置信度过低时：
•	FacePresence(k) = 0
•	Face\_raw(k) 维持上一有效值或全片中位数，保持曲线连续；后续通过 mask 强制归零

最终：
[
FaceAffect(k) = \sigma(\text{zscore}(\text{Face\_raw}(k))), \quad FacePresence(k) \in \{0,1\}
]
并在 Arousal 中使用：
[
A_{face}(k) = FacePresence(k) \cdot FaceAffect(k)
]
FacePresence 以单独曲线形式写入数据库，供 Arousal、UI、导出使用；若 4.5.2 判定整片无主角轨迹，则仅写 FacePresence=0，Arousal 中删除 face 维度。

4.5.4 范围说明
ArcScope 显示的是电影的视听节奏，而非心理学。对白语义、群体心理仍由导演/分析者自行解读；系统只提供：
•	主角情绪轨迹（FaceAffect）
•	整体张力（Pace/Sound/Motion/Info）
吵架戏/群戏的语义差异通过声音强度、语速、Motion、景别等指标体现，而非靠平均每张脸的情绪。
若“无可靠主角轨迹”，ArcScope 直接关闭 FaceAffect 曲线与 FacePresence mask（全 0），UI/报告均提示 Arousal 仅由 Pace/Sound/Info 驱动，确保系统稳定而非勉强输出噪声。
 
4.6 情绪 / 唤醒曲线 Arousal(t)
4.6.1 含义
整合节奏、声音、运动、表演等多模态信息，给出每秒的“兴奋度/紧绷度” proxy。
4.6.2 输入特征
已得到的曲线：
•	Pace(t_k)
•	SoundEnergy(t_k)
•	ColorMood(t_k)
•	InfoDensity(t_k)
•	FaceAffect(t_k)
另外加入原始运动强度 (M(k)) 的变化率：
[
\Delta M(k) =
\begin{cases}
0, & k=0 \
|M(k) - M(k-1)|, & k \ge 1
\end{cases}
]
所有序列直接使用其稳健 z-score 形式 (Z_{\bullet})。
4.6.3 曲线定义：多模态隐变量
将上述特征堆叠为矩阵：
[
\mathbf{Y}(t) = [Z_{Pace}(t), Z_{Sound}(t), Z_{Color}(t), Z_{Info}(t), Z_{Face}(t), Z_{\Delta M}(t)]
]
ArcScope 对 \mathbf{Y}(t) 做 PCA，取第一主成分 (w_{arousal}) 作为该片的唤醒主方向：
[
\text{Arousal_raw}(t) = \mathbf{Y}(t) \cdot w_{arousal}
]
为确保“高能量”方向一致，可用 Pace/Sound 高段的符号作为引用，再稳健 z-score + sigmoid：
[
Arousal(t_k) \in (0,1)
]
其中面部维度使用 A_{face}(t) = FacePresence(t) · \widetilde{FaceAffect}(t)：mask 为 0 时 face 贡献强制为 0；mask 为 1 时正常参与 PCA。若整片没有可靠主角轨迹，则 FacePresence 全为 0，Arousal 永远落在 Pace/Sound/Info/Motion 的子空间。
4.6.4 高级版（预留）
未来如需替换，可在此 PCA 基线之上引入更复杂模型；1.0 固定采用 PCA PC1。

4.6.5 情绪累积进度（跨片对齐用）
为了在不同片长之间对齐“情绪旅程”，ArcScope 以 Arousal(t) 定义情绪进度：
1.	累计唤醒量：
[
S(t) = \int_{0}^{t} A(\tau)\, d\tau
]
2.	情绪进度（归一化到 [0,1]）：
[
u(t) = \frac{S(t)}{S(T)} = \frac{\int_{0}^{t} A(\tau)\, d\tau}{\int_{0}^{T} A(\tau)\, d\tau}
]
离散情况下，以每秒累加 A(k) 并除以全片累积和即可；u=0 表示情绪旅程刚开始，u=1 表示完成。该 u(t) 仅用于跨影片对齐，单片浏览仍使用线性时间 t 与归一化时间 \tau = t / T。

4.7 多曲线耦合的自适应建模
ArcScope 不根据“类型”写死规则，而是让多模态特征通过 PCA 自动竞争：
•	对 Pace/Sound/Color/Info/Face 矩阵做 PCA，第一主成分即该片在该子空间的主方向
•	若影片由 Motion 主导，PC1 自然包含更多运动；剧情片会让 Info loading 更高
•	参考曲线 overlay 亦可将经典影片的 PCA 方向存档，在共享主轴中比较
整套系统保持线性、可解释、无人工权重。
 
5. 多尺度结构分段（Scenes / Sequences / Acts）
5.1 连续曲线与离散块的统一
ArcScope 的全部输入都是连续曲线 f(t)，而多尺度结构（场景、序列、幕）是一系列时间段标签。二者不矛盾：
•	曲线提供“每秒发生了什么”的定量信号
•	分段结构是在这些信号上“切块”，标注某段时间属于哪个 Scene / Sequence / Act
因此无需新增曲线，只需在现有时间轴上自动生成分段色带即可。

5.2 Scene 级别：节奏-情绪耦合的变点检测
•	基础信号：x(t) = [Pace, Sound, Color, Info]
•	情绪动力：使用 3.3.2 的 A'(t)、A''(t)、D_{PA}(t)
•	综合变点指标：
[\
S(t) = w_1 ∥x(t) - x(t-1)∥_2 + w_2 |A'(t)| + w_3 |A''(t)| + w_4 |D_{PA}(t)|
]
其中 w = [w_1…w_4] 由对 [∥x差值∥, |A'|, |A''|, |D_{PA}|] 做 PCA 得到（PC1 loadings）。
•	窗口与间距：使用 3.3 的 W，Scene 边界需间隔 ≥ W/2
•	变点判定：S(t) > Q_S(0.85) 的时刻标记为 scene change（分位数自适应）

5.3 Sequence / Act 级别：联合聚合
•	Scene 指纹：
[\
v_i = [\text{mean}(x), \text{mean}(A), \text{mean}(A'), \text{mean}(A''), \text{mean}(D_{PA})]
]
•	仅合并相邻 Scene，距离 d(i,j) = ∥v_i - v_j∥_2
•	Sequence 目标：K_seq = min(12, max(4, ⌊N_scene/3⌋))，greedy 合并最小代价对直至剩 K_seq 段
•	Act 目标：总时长 <60 min →2；60–150 min →3；≥150 min →4；在 Sequence 级重复相邻合并
•	场景色彩指纹继续记录 Hue 直方图、ColorWarmth/Energy/Harmony 均值，用于 Poster/metadata

5.4 Scene 类型自动标签（情绪-结构语义）
根据 Scene 指纹自动打标签（阈值均为片内分位数）：
1.	Low-Energy Basin：mean(A) < Q_A(0.3)，mean(Pace) < Q_P(0.3)
2.	Build-up：mean(A') > Q_{A'}(0.7)，mean(D_{PA}) < Q_{D_{PA}}(0.3)
3.	High-Energy Peak：mean(A) > Q_A(0.75)
4.	Emotional Reversal：|mean(A')| > Q_{|A'|}(0.8) 且 valence 跨越中点
5.	Misaligned Tension：mean(D_{PA}) > Q_{D_{PA}}(0.8)
6.	Release：mean(A') < Q_{A'}(0.3)，且 Pace/Sound 下行
这些标签在 UI/报告中显示，帮助理解段落语义。

5.5 可视化：曲线时间线下的多层色带
•	在 CurvesTimelineView 下方叠加多条色条轨道：Shot 网格 / Scene 色块 / Sequence 色块 / Act 色块
•	Scene 轨道色块较密，Sequence 稍粗，Act 最粗，形成“嵌套”视觉
•	点击任意色块可联动曲线高亮对应区间，也可与 Diagnostics Heatmap 共用交互
•	导出 Poster 时，将 Scene/Sequence/Act 轨道一并绘制，帮助用户一眼看到结构

5.6 情绪进度对齐 + 实现
•	Scene/Sequence 切分在执行前可选地对所有曲线进行情绪进度重参数化 u(t)（4.6.5），即在等情绪步长的坐标上检测变点，实现“情绪旅程”对齐的结构切片
•	FilmEngine 在生成曲线后即刻运行上述 S(t)+greedy 聚合流程，输出 scene/sequence/act 三张 segment 表
•	SQLite 可复用 segments 表或新增 scene_segments / act_segments 以存储 start/end/level、scene_type、色彩/情绪指纹
•	SwiftUI 读取这些分段，与 Diagnostics 共享时间轴
•	FilmEngine 在生成曲线后即刻运行上述 D(t)+greedy 聚合流程，输出 scene/sequence/act 三张 segment 表
•	SQLite 可复用 segments 表或新增 scene_segments / act_segments 以存储 start/end/level
•	SwiftUI 读取这些分段，与 Diagnostics 共享时间轴

5.6 Metadata Overlays（展示修辞，不入模型）
ArcScope 在主时间线上提供以下可选轨道，默认隐藏，开启后仅作语义注释：
•	Audio：BPM(t)、TempoConfidence(t)、Key(t)/KeyChangeDensity(t)、HarmonicTension(t)、StereoWidth(t)、ReverbAmount(t)、SilenceScore(t)（曲线）
•	Photography/Color：ShotScale 色带（ELS→ECU）、ShotScaleNumeric 曲线、CameraMotionType 色带（推/拉/摇/移/跟）、Camera vs Object Motion 对比曲线、Hue Track / Color Harmony (ΔE) 条带
•	Face/Dialogue：LookAtCamera(t)、FrontalScore(t)、DowncastScore(t)、IsSpeaking(t)、Dialogue/Music/Effects Dominance 0/1 条带
•	Info：ExpositionCurve(t)（解释性比重）、Verbal/Visual/Event Load 子曲线
这些 overlay 永不进入 PCA/Arousal/Diagnostics，只提供阅读上下文。

6. 自动诊断引擎设计
6.1 数据结构
enum class IssueType {
    LowActivity,
    OvercutFlat,
    AudioVisualMisalign,
    InfoKillsEmotion,
    HighTensionLowPace,
    SoundOverdriveFlatEmotion,
    SoundColorConflict,
    ColorEmotionMismatchWarm,
    ColorEmotionMismatchCold,
    CognitiveOverload,
    FacialAffectSuppressed,
    NonCharacterDrivenPeak,
    SilentHighTension,
    ModeCounterpoint
};

struct IssueSegment {
    IssueType type;
    double    start_sec;
    double    end_sec;
    float     severity;  // 0–1，基于偏离程度
};
6.2 基本思路
•	所有曲线均在 ([0,1])，先计算关键分位数：Q_X(p) = 第 p 分位点
•	持续时间不写死秒数，而是根据每条曲线的“稳定段”中位长度 (L_X) 决定，示例：
[\
T_{min}^{(P)} = \operatorname{median}\big\{ \text{length of segments where } P(t) > Q_P(0.5) \big\}
]
•	异常程度统一按“偏离分位点的归一面积”计算，而不是绝对差
•	每个 issue 都只依赖影片自身分布，从而对快剪/慢片都成立
6.3 规则示例
记：
•	P(t) = Pace(t)
•	S(t) = SoundEnergy(t)
•	C(t) = ColorMood(t)
•	I(t) = InfoDensity(t)
•	A(t) = Arousal(t)
•	F(t) = FaceAffect(t)（仅主角在场时定义，否则视为缺失/0）

共用量（全系统统一）：
•	High_X = Q_X(0.75)、VeryHigh_X = Q_X(0.8)
•	Low_X = Q_X(0.25)、VeryLow_X = Q_X(0.3)
•	T_{min}^{(X)} = median stable length of X（按 6.2 公式）

6.3.1 LowActivity（低活跃盆地）
条件：
[
P(t) < Low_P,\ S(t) < Low_S,\ A(t) < Low_A,\ I(t) < Low_I
]
持续时间 ≥ T_{min}^{(P)}（若影片整体平缓，T_{min} 自动较长；快片则较短）。
严重度：
[
severity = \frac{1}{t_b - t_a} \int_{t_a}^{t_b} \Big( \frac{Low_P - P(t)}{Low_P} + \frac{Low_S - S(t)}{Low_S} + \frac{Low_A - A(t)}{Low_A} + \frac{Low_I - I(t)}{Low_I} \Big) dt
]

6.3.2 OvercutFlat（节奏过猛但情绪空转）
条件：
[
P(t) > High_P,\quad A(t) < Low_A
]
持续时间 ≥ T_{min}^{(P)} / 2（因为快节奏段天然较短，T_{min} 自动缩短）。
严重度以 Pace 超过 High_P 的幅度和 Arousal 低于 Low_A 的幅度加权平均。

6.3.3 AudioVisualMisalign（视听错位）
在计算前，先估计全局 AV 延迟：
•	取 Pace/SoundEnergy 的 z-score 序列 Z_P(k)、Z_S(k)
•	在 ℓ ∈ {-3,…,3} 秒内搜索最大化 R(ℓ) = \sum_k Z_P(k)·Z_S(k+ℓ) 的偏移 ℓ*
•	后续所有比较均使用校正后的 Z_S(k+ℓ*)

定义声音与画面节奏差（使用中心对齐的 t_k = k + 0.5 秒）：
[
D_{PS}(t) = Z_P(t) - Z_S(t)
]
其中 Z 为对应曲线的稳健 z-score。找出 |D_{PS}(t)| > Q_{|D|}(0.8) 且持续 ≥ T_{min}^{(S)} 的区段；
正差 → 画面快但声音弱，负差 → 声音炸但画面平。

6.3.4 InfoKillsEmotion（解释压死情绪）
条件：I(t) > Q_I(0.8)、A(t) < Q_A(0.3)、Exposition(t) > Q_{Expo}(0.6)
持续：duration ≥ median length of Info 高段
severity = \int (I-High_I)_+ dt · \int (Low_A - A)_+ dt / L

6.3.5 HighTensionLowPace（静态高张力）
条件：A(t) > Q_A(0.85) 且 P(t) < Q_P(0.25)
持续：duration ≥ 0.6·T_{min}(P)
severity = mean\Big(\frac{A(t)-Q_A(0.85)}{Q_A(0.85)}\Big)

6.3.6 SoundOverdriveFlatEmotion（声音驱动失败）
条件：SoundEnergy(t) > Q_S(0.75) 且 A(t) < Q_A(0.25)
持续：duration ≥ T_{min}(S)
severity = mean(Z_S(t) + Z_{LowA}(t))

6.3.7 SoundColorConflict（声音与色彩冲突）
定义 ColdDark(t) = Z_{Gray}(t) + Z_{CoolHue}(t)
条件：SoundEnergy(t) > Q_S(0.75) 且 ColdDark(t) > Q_{Cold}(0.75)
持续：duration ≥ min_dur(S, ColdDark)
severity = mean(Z_S(t) + Z_{Cold}(t))

6.3.8 ColorEmotionMismatch（色彩与情绪背离）
使用 AlignNorm(t)（4.3.4）。
条件：AlignNorm(t) > Q_{Align}(0.8)，或 (A(t) > Q_A(0.7) 且 Brightness(t) < Q_B(0.3))，或 (A(t) < Q_A(0.3) 且 ColdDark(t) > Q_{Cold}(0.7))
持续：duration ≥ T_{min}(A)
severity = mean(AlignNorm(t))。可区分暖亮/冷暗两类 IssueType。

6.3.9 CognitiveOverload（认知超载）
条件：I(t) > Q_I(0.75) 且 P(t) > Q_P(0.75)
持续：duration ≥ min_dur(I,P)
severity = mean(Z_I(t) + Z_P(t))

6.3.10 FacialAffectSuppressed（表演被压制）
条件：FaceAffect(t) > Q_F(0.75) 且 A(t) < Q_A(0.25)
持续：duration ≥ T_{min}(Face)
severity = mean(Z_{Face}(t) - Z_A(t))

6.3.11 NonCharacterDrivenPeak（非角色情绪高潮）
条件：A(t) > Q_A(0.85) 且 FaceAffect(t) < Q_F(0.25)
持续：duration ≥ 0.5·T_{min}(A)
severity = mean(Z_A(t) - Z_{Face}(t))

6.3.12 SilentHighTension（静默高张力）
条件：SoundEnergy(t) < Q_S(0.3) 且 Tension(t) > Q_T(0.7)
持续：duration ≥ T_{min}(Tension)
severity = mean(Z_Tension(t) - Z_S(t))

6.3.13 ModeCounterpoint（多模态反向）
定义 sign_X(t) = sign(Z_X(t))。
若 Pace/Sound/Color 至少两者与 Arousal 同时符号相反，且持续 ≥ 0.5·W，则判定。
severity = 冲突时刻占窗口比例。

6.4 Heatmap 颜色编码
CurvesTimelineView 底部的 Diagnostics Heatmap 使用固定色相区分 Issue 类型，深浅表示 severity：
•	LowActivity：淡蓝→深蓝
•	OvercutFlat / HighTensionLowPace：淡橙→深橙
•	AudioVisualMisalign / SoundOverdrive：淡紫→深紫
•	InfoKillsEmotion / CognitiveOverload：淡红→深红
•	ColorEmotionMismatch / SoundColorConflict / ModeCounterpoint：青色或青紫系（按启用 issue 扩展）
交互规则：悬停/点击时显示颜色、icon（“…” 代表信息、🔊 代表声音等）、文本说明，保持一致性、不混用色相。
7. FilmEngine 设计（C++17）
7.1 核心类型
enum class CurveType {
    Pace,
    Sound,
    ColorMood,
    InfoDensity,
    Arousal,
    FaceAffect
};

struct Curve {
    CurveType           type;
    double              fps;    // 通常 1.0
    std::vector<float>  values; // 归一化到 0~1
};

struct FilmCurves {
    Curve pace;
    Curve sound;
    Curve color_mood;
    Curve info;
    Curve arousal;
    Curve face_affect;    // 可能为空（无主角轨迹时不写入）
    Curve face_presence;  // 0/1 掩码，长度与其他曲线一致
};

struct FilmAnalysisResult {
    FilmCurves                curves;
    std::vector<IssueSegment> issues;
};
7.2 模块划分
1.	MediaReader
o	负责基本时间轴信息（影片时长）
o	实际帧/音频读取由 Bridge 使用 AVFoundation 完成，C++ 层主要使用 Bridge 提供的特征数据；如果需要，也可预留 FFmpeg 后备实现。
2.	features::SoundFeatures
o	输入：音频 PCM 分块
o	输出：每秒特征
    	RMS 能量 E_rms
    	低频 / 高频能量 E_LF / E_HF
    	节奏：Spectral Flux、局部 BPM、TempoConfidence、RhythmEnergy
    	和声：Chroma 向量、Key(k)、HarmonicTension
    	可选：speech probability / VAD（供 InfoFeatures 复用）
3.	features::MotionFeatures & features::ShotFeatures
o	输入：每秒光流平均幅度（来自 Vision）
o	输出：运动强度 (M(k))、镜头边界、local ASL、cut density
4.	features::ColorFeatures
o	输入：每秒采样帧（RGB）
o	输出：(L_{mean}, S_{mean}, Warmth, Grayness)
5.	features::InfoFeatures
o	输入：字幕列表 + ASR/VAD 结果 + SoundEvents 检测 + 视觉熵统计 + Face 模块输出
o	输出：每秒 Words、SpeechDuty、MusicChange、SoundEvents、VisualEntropy、FaceChange，以及 ExpoRatio（metadata 用）
6.	features::FaceFeatures
o	输入：Bridge 输出的人脸检测 + embedding + CoreML 表情结果
o	流程：
    	多目标跟踪生成 track_i，累计 T_i、S_i、Q_i，按 Score_i = T_i·S_i·Q_i 选出主角轨迹
    	在每秒时间轴上插值出主角 presence/面积/表情/valence/arousal，并标记主角缺席区间
o	输出：
    	DominantFaceCurve（主角轨迹 + 置信度 + 景别权重）
    	FacePresenceMask（布尔/权重，用于 FaceAffect/Arousal gating）
    	无可靠主角 → 输出空 DominantFaceCurve、全 0 mask，并向 FilmAnalyzer 报告“face disabled”
7.	CurveEngine
o	组合上述特征，按第 4 章定义生成 5–6 条曲线。
8.	DiagnosticsEngine
o	扫描曲线，输出问题区段列表。
9.	FilmAnalyzer
class FilmAnalyzer {
public:
    FilmAnalyzer(const std::string& path, double fps = 1.0);

    FilmAnalysisResult run(ProgressCallback cb = nullptr);
private:
    std::string path_;
    double      fps_;
};
 
8. 数据与持久化（SQLite）
8.1 Schema 设计原则

ArcScope 采用"逻辑 ID + 外部对齐 key"的双重标识策略：

1.	logical_id：稳定的用户可见标识符（如 blade-runner-1982），作为主键
2.	外部对齐 key：可选的 imdb_id、tmdb_id、letterboxd_id 等，用于跨平台数据对齐
3.	file_path：实际文件路径（可空，本地分析时填写）
4.	source：数据来源（local / imdb / omdb / tmdb 等）

关键设计要点：
•	logical_id 必须唯一，作为影片的主键（UNIQUE）
•	imdb_id 非唯一（允许同一 IMDb ID 对应多个版本/剪辑）但建议加 INDEX
•	updated_at 字段必须存在，用于按最近分析排序（ORDER BY updated_at DESC）
•	UPSERT 策略：按 logical_id 冲突更新，而非 id

8.1.1 films 表（核心）
CREATE TABLE IF NOT EXISTS films (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  logical_id    TEXT NOT NULL UNIQUE,              -- 稳定的逻辑标识符（如 blade-runner-1982）
  title         TEXT NOT NULL,
  year          INTEGER,
  director      TEXT,
  country       TEXT,
  runtime_sec   REAL,
  file_path     TEXT,                              -- 实际文件路径（本地分析时填写）
  imdb_id       TEXT,                              -- IMDb ID（如 tt0083658），可选
  tmdb_id       INTEGER,                           -- TMDb ID（可选）
  source        TEXT,                              -- 数据来源：local / imdb / omdb / tmdb 等
  tags          TEXT,
  my_rating     REAL,
  ext_rating    REAL,
  created_at    TEXT DEFAULT CURRENT_TIMESTAMP,
  updated_at    TEXT DEFAULT CURRENT_TIMESTAMP     -- ⚠️ 必须存在（ORDER BY updated_at）
);

CREATE INDEX IF NOT EXISTS idx_films_imdb ON films(imdb_id);      -- 支持 IMDb 查询
CREATE INDEX IF NOT EXISTS idx_films_logical ON films(logical_id); -- 冗余但加速查询
CREATE INDEX IF NOT EXISTS idx_films_updated ON films(updated_at); -- 支持"最近分析"排序

-- ⚠️ TRIGGER：自动更新 updated_at
CREATE TRIGGER IF NOT EXISTS update_films_timestamp
  AFTER UPDATE ON films FOR EACH ROW
BEGIN
  UPDATE films SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
END;

8.1.2 IMDb ID 集成策略

能否直接用 IMDb ID 当 film_id？

推荐：支持但不强制。使用"逻辑 ID 主导 + IMDb ID 辅助"的策略。

实施建议：
1.	保留 logical_id（如 blade-runner-1982）作为主键
2.	新增 imdb_id 字段（如 tt0083658）作为外部对齐 key
3.	UI 导入时：
    •	用户可输入/匹配 IMDb ID → 系统自动填充 title/year/director
    •	也可纯本地导入（无 IMDb ID）→ 系统从文件名/元数据生成 logical_id
4.	多外部 ID 并存：未来可扩展 tmdb_id、letterboxd_id、douban_id 等

优势：
•	避免对 IMDb 的硬依赖（短片/实验片/私人剪辑版本可能没有 IMDb）
•	支持同名/重名影片（通过 logical_id 区分）
•	多平台数据对齐（Letterboxd / TMDb / Douban）

OMDb API 集成示例（推荐）：
// 用户输入片名/年份 → 调用 OMDb API
// GET http://www.omdbapi.com/?apikey=YOUR_KEY&t=Blade+Runner&y=1982
// 返回 JSON 包含 imdbID: "tt0083658"
// 存储到 films.imdb_id

注意：OMDb 免费层有速率限制（1000 req/day），生产环境需付费或缓存。

8.1.3 其他表

CREATE TABLE IF NOT EXISTS curves (
  id         INTEGER PRIMARY KEY AUTOINCREMENT,
  film_id    INTEGER NOT NULL,
  curve_type TEXT NOT NULL,   -- 'pace','sound','color','info','arousal','face'
  fps        REAL NOT NULL,
  length     INTEGER NOT NULL,
  data_blob  BLOB NOT NULL,
  FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_curves_film ON curves(film_id);

约定：FaceAffect 与 FacePresence 分别占据一行（curve_type = 'face_affect', 'face_presence'），二者共享相同长度；FacePresence 仅存 {0,1} 掩码。若整片无主角轨迹，仅写入 face_presence（全 0），省略 face_affect 行。

CREATE TABLE IF NOT EXISTS metrics (
  id          INTEGER PRIMARY KEY AUTOINCREMENT,
  film_id     INTEGER NOT NULL,
  name        TEXT NOT NULL,
  value       REAL NOT NULL,
  FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS segments (
  id          INTEGER PRIMARY KEY AUTOINCREMENT,
  film_id     INTEGER NOT NULL,
  type        TEXT NOT NULL,   -- 'low_activity','overcut_flat','av_misalign','info_kill_emotion'
  start_sec   REAL NOT NULL,
  end_sec     REAL NOT NULL,
  severity    REAL NOT NULL,
  FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS references_library (
  id          INTEGER PRIMARY KEY AUTOINCREMENT,
  film_id     INTEGER NOT NULL,
  note        TEXT,
  FOREIGN KEY(film_id) REFERENCES films(id) ON DELETE CASCADE
);
8.2 关键实现 Bug 与修复清单（⚠️ 必须在 1.0 前修复）

8.2.1 BUG-001：films 表缺少 updated_at 列
位置：core/src/SQLiteStore.cpp:74-89
影响：arcscope_bridge.mm:234 中 ORDER BY updated_at DESC 会导致 SQLite 错误

当前 schema（错误）：
CREATE TABLE IF NOT EXISTS films (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  ...
  created_at TEXT DEFAULT CURRENT_TIMESTAMP  -- ❌ 缺少 updated_at
);

修复方案：
1. 在 schema 中添加 updated_at 列
2. 添加 TRIGGER 自动更新该字段（见 8.1.1）
3. 为已部署数据库添加 migration：
   ALTER TABLE films ADD COLUMN updated_at TEXT DEFAULT CURRENT_TIMESTAMP;

8.2.2 BUG-002：SQLiteStore::save_analysis 的 UPSERT 逻辑错误
位置：core/src/SQLiteStore.cpp:192-208
影响：重复导入影片时会创建新记录，而非更新现有记录

当前代码（错误）：
auto stmt = prepare_or_throw(db_,
    "INSERT INTO films (title, year, director, runtime_sec, file_path) VALUES (?,?,?,?,?) "
    "ON CONFLICT(id) DO UPDATE SET ..."   // ❌ id 从未提供，UPSERT 永远不会触发
    "RETURNING id;");

sqlite3_bind_text(stmt, 5, result.filmId.c_str(), ...);  // ❌ filmId 绑定到 file_path，语义错位

问题分析：
1. ON CONFLICT(id) 需要在 INSERT 时显式提供 id，但当前代码从未提供
2. result.filmId 是逻辑 ID（如 blade-runner-1982），应绑定到 logical_id，而非 file_path
3. file_path 应该存储实际文件路径（如 /Users/foo/Videos/blade_runner.mp4）

修复方案：
auto stmt = prepare_or_throw(db_,
    "INSERT INTO films (logical_id, title, year, director, runtime_sec, file_path, imdb_id, source, updated_at) "
    "VALUES (?,?,?,?,?,?,?,?, CURRENT_TIMESTAMP) "
    "ON CONFLICT(logical_id) DO UPDATE SET "  // ✅ 按 logical_id 冲突
    "title=excluded.title, year=excluded.year, director=excluded.director, "
    "runtime_sec=excluded.runtime_sec, file_path=excluded.file_path, "
    "imdb_id=excluded.imdb_id, source=excluded.source, updated_at=CURRENT_TIMESTAMP "
    "RETURNING id;");

sqlite3_bind_text(stmt, 1, result.filmId.c_str(), ...);       // logical_id
sqlite3_bind_text(stmt, 2, result.filmTitle.c_str(), ...);    // title
sqlite3_bind_int(stmt, 3, result.year);                       // year
sqlite3_bind_text(stmt, 4, result.director.c_str(), ...);     // director
sqlite3_bind_double(stmt, 5, result.durationSeconds);         // runtime_sec
sqlite3_bind_text(stmt, 6, result.filePath.c_str(), ...);     // file_path (实际路径)
sqlite3_bind_text(stmt, 7, result.imdbId.c_str(), ...);       // imdb_id (可空)
sqlite3_bind_text(stmt, 8, result.source.c_str(), ...);       // source

8.2.3 BUG-003：Bridge 层时间戳不一致（Contract 违反）
位置：bridge/VideoFeatureExtractor.mm:427 vs bridge/arcscope_bridge.mm:165-174, 479
影响：Bridge 输出的 time_seconds = i，但 Core 和 UI 期望 i+0.5（中心对齐）

Contract 定义（arcscope_bridge.mm:165）：
// TIME INVARIANT CONTRACT (§ I.R3):
// samples[i].time_seconds MUST equal i + 0.5 (center-aligned 1Hz)

当前实现（错误）：
// VideoFeatureExtractor.mm:427
sample.time_seconds = static_cast<double>(i);  // ❌ 左端点 (i)，违反 contract

// arcscope_bridge.mm:479（读取曲线）
curve.times[i] = static_cast<float>(i) + 0.5f; // ✅ 中心对齐 (i+0.5)

不一致后果：
1. UI 曲线显示与 FilmEngine 内部计算的 t_center 错位 0.5 秒
2. 跨影片参考曲线对齐时会出现系统性偏移
3. Scene/Sequence 边界检测的时间点不准确

修复方案：
// VideoFeatureExtractor.mm
for (size_t i = 0; i < sampleCount; ++i) {
    ArcScopeBridgeSample sample;
    sample.time_seconds = static_cast<double>(i) + 0.5;  // ✅ 强制中心对齐
    ...
}

8.2.4 BUG-004：FilmAnalysisResult 缺少必要字段
位置：core/include/arcscope/FilmTypes.h
影响：save_analysis 无法正确写入 logical_id、file_path、imdb_id、source

当前结构（不完整）：
struct FilmAnalysisResult {
    std::string filmId;      // 实际是 logical_id
    std::string filmTitle;
    int year{0};
    std::string director;
    double durationSeconds{0.0};
    // ❌ 缺少 filePath、imdbId、source
    ...
};

修复方案：
struct FilmAnalysisResult {
    std::string filmId;         // logical_id（如 blade-runner-1982）
    std::string filmTitle;      // 显示标题
    std::string filePath;       // ✅ 实际文件路径（如 /Users/foo/Videos/...）
    std::string imdbId;         // ✅ IMDb ID（如 tt0083658），可空
    std::string source;         // ✅ 数据来源（local / imdb / omdb）
    int year{0};
    std::string director;
    double durationSeconds{0.0};
    ...
};

8.2.5 修复优先级

P0（阻塞 1.0）：
1. BUG-002：修复 UPSERT 逻辑（否则重复导入会失败）
2. BUG-003：修复时间戳对齐（否则所有时间相关功能都会错位）
3. BUG-001：添加 updated_at 列（否则列表排序会崩溃）

P1（影响体验）：
4. BUG-004：补全 FilmAnalysisResult 字段（否则无法存储完整元数据）

8.3 封装接口（简要）
namespace arcscope::db {

class Connection {
public:
    explicit Connection(const std::string& path);
    ~Connection();
    sqlite3* handle() const;
};

struct FilmRecord {
    int         id;
    std::string title;
    int         year;
    std::string director;
    double      runtime_sec;
};

int  insert_film(Connection&, const FilmRecord&);
void insert_curves(Connection&, int film_id, const FilmCurves&);
void insert_metrics(Connection&, int film_id, const FilmCurves&);
void insert_segments(Connection&, int film_id,
                     const std::vector<IssueSegment>&);

std::vector<FilmRecord> list_films(Connection&);

}
 
9. Bridge 层设计（Objective-C++）
9.1 C 接口（供 Swift 调用）
typedef void (*ArcScopeProgressCallback)(double progress);

int arcscope_analyze_film(
    const char* input_path,
    const char* db_path,
    const char* title,
    int         year,
    const char* director,
    ArcScopeProgressCallback cb);

typedef struct {
    int         id;
    const char* title;
    int         year;
    const char* director;
    double      runtime_sec;
} ArcScopeFilmSummary;

int arcscope_list_films(const char* db_path,
                        ArcScopeFilmSummary* out_array,
                        int max_count);
9.2 职责
•	用 AVFoundation 打开视频，获取时长、基础元数据。
•	调用 Vision：
o	VNGenerateOpticalFlowRequest 生成光流，按秒统计 flow_mag。
o	人脸检测 + 人脸 ROI 提取。
•	调用 CoreML 表情模型：
o	得到每秒 face presence / intensity / valence / arousal。
•	汇总为 features::MotionFeatures 与 features::FaceFeatures 传给 C++。
•	调用 FilmAnalyzer::run，传入进度回调。
•	调用 db 封装，将结果写入 SQLite。
 
10. macOS App（SwiftUI）设计
ArcScope macOS App 采用 SwiftUI + ObservableObject 架构，整体风格为 Cinematic Dark Mode。模块包含 LibraryView、FilmDetailView、ReferencePicker、SettingsView、PosterExportView 等。所有数据通过 SQLite→ViewModel→SwiftUI 的单向流。

10.1 LibraryView
职责：展示影片库、触发新分析。结构：NavigationSplitView（Sidebar + detail）。Sidebar 列表使用 FilmCard（暗色卡片，含缩略图/标题/年份/时长），toolbar 提供“Analyze New Film”按钮；detail 根据 selectedFilm 显示 FilmDetailView 或提示文本。LibraryViewModel 管理 films、selectedFilm、importFilm()（调用 arcscope_analyze_film）。

10.2 FilmDetailView
职责：展示单片所有曲线、诊断、结构色带及参考 overlay。VStack(spacing:0)：FilmHeaderView（标题/年份/ReferencePicker/导出按钮）→ CurvesTimelineView（曲线主画布）→ DiagnosticsHeatmap → StructureBands。背景 Color.black.opacity(0.92)，onAppear 加载曲线。

10.3 CurvesTimelineView
多层 ZStack：GradientBackdrop + GridOverlay + CurveLayers（Pace/Sound/Color/Info/Arousal/Face，Catmull-Rom Path、渐变填充、glow）+ ReferenceOverlayLayer（使用情绪进度 u(t) 重采样的幽灵曲线）+ HUDCursor（显示当前时间点各曲线数值）+ GestureLayer（DragGesture 更新 cursorX）。

10.4 DiagnosticsHeatmap & StructureBands
DiagnosticsHeatmap 通过 GeometryReader + ForEach(issue) 绘制条带，颜色遵循 6.4。StructureBands 在 Scene/Sequence/Act 三层绘制矩形，显示 start/end、scene_type label、色彩指纹；背景半透明黑。

10.5 ViewModel 逻辑
• LibraryViewModel：films、selectedFilm、importFilm
• FilmCurvesViewModel：加载曲线/诊断/结构/metadata，提供 value(at:)、x(forTime:)、reference 重采样（基于 u(t)）等 API，并管理曲线/metadata overlay toggles。

10.6 Reference / Settings / Poster Export
ReferencePickerView 供选择参考影片；SettingsView（Form + AppStorage）配置分析参数；PosterExportView 展示 Poster 布局（曲线缩略、Keyframe、色彩 fingerprint），使用 ImageRenderer 导出 PDF/PNG。

10.7 Metadata Overlays
ShotScale、CameraMotion、FacePresence、ExpositionCurve、Verbal/Visual/Event Load、Hue Track、ColorHarmony 等以独立轨道展示（默认隐藏，toggle 控制），渲染方式与 StructureBands 类似。

10.8 SwiftUI 技巧
• 曲线绘制封装 Path + Catmull-Rom，并通过 .shadow/.fill 增加荧光效果
• HUDCursor 使用 GeometryReader + DragGesture + .ultraThinMaterial 背景
• Heatmap/Bands 等全用 ForEach + offset，便于扩展
• 数据以 ObservableObject/StateObject 注入视图；参考 overlay 通过 u(t) 重采样
• Poster 导出采用 ImageRenderer；SettingsView 使用 Form + AppStorage

10.9 UI 风格：Cinematic Dark Mode
•	默认采用沉浸式暗色模式，背景使用深灰蓝（例如 Color(red: 0.1, green: 0.1, blue: 0.12)），避免纯黑造成眩光
•	全局界面引用 NSVisualEffectView / .background(.ultraThinMaterial)，让影片色块透入侧栏与浮动面板，制造影院般的层次
•	曲线颜色遵循“荧光”调性，配合深背景产生发光感；辅助网格线极淡，仅指引而不抢戏
•	数字、坐标轴、HUD 使用 SF Mono，正文和按钮沿用 SF Pro，营造“精密仪器”气质
•	若检测不到主角轨迹，FaceAffect/FacePresence 的开关与图例默认隐藏，并在信息栏提示 “No stable dominant face; arousal driven by pace+sound”

10.4 Curve Timeline 设计（灵魂界面）
•	曲线绘制采用层叠面积 / ridge plot 变体：
    o	Pace：亮橙线 + 50%→0% 渐变填充，如火焰
    o	Sound：可绘制为示波器式实心波形或荧光线
    o	Color/Info：用半透明色带区分，避免毛线球
    o	Arousal：加粗红线，高值处叠加 .shadow(color: .red, radius: 5) 的辉光
•	交互式游标：拖动时在当前秒显示 HUD，列出 Pace/Sound/... 的即时值，用户无需回看 Y 轴
•	问题区段：时间轴底部绘制热力图带，正常透明，按类型使用固定色相（蓝/橙/紫/红）加深浅表示 severity，类似现代股票图的 warning stripe
•	多尺度结构：Scene/Sequence/Act 色带与 Heatmap 同轴显示，可展开/收起，方便把曲线与自动切块对齐
•	Cinematography Metadata：可选轨道（ShotScale/CameraMotion/Speaking vs Listening），默认隐藏，启用后以色条或 0/1 条带展示，仅作语义注释
•	FacePresence mask：面部曲线保持连续；mask=0 的时间段用虚线/低透明显示，并可在最底部添加 0/1 细轨辅助调试
•	HUD/Legend 在 mask=0 时提示 “Main character off-screen / low confidence”，便于用户理解 face 贡献为何被禁用
•	Heatmap 配色遵循文档 6.4：蓝=LowActivity、橙=OvercutFlat、紫=AV Misalign、红=Info，深浅代表 severity，鼠标 hover 展示 icon+文本
•	参考金融/音频仪表的渐变与描边，保持高级而克制

10.5 导出“Souvenir”报告
•	导出模式 A：Film Fingerprint Poster
    o	顶部：电影海报 + 片名/导演/年份/时长，套用 Criterion 式极简排版
    o	中部：5 条核心曲线分轨排布，仿心电图拉平（不重叠）
    o	底部：在 Pace/Arousal 峰值自动抓取关键帧缩略图，标注时间戳
•	导出模式 B：标准分析报告（PDF/Markdown）沿袭 UI 色彩，方便论文/社媒分享
•	“Souvenir” 输出即免费广告，图形要有装饰性，可直接打印装裱

10.6 参考影片的“幽灵模式”
•	当叠加参考曲线（如《教父》）时，以虚线或低透明灰色填充呈现；当前影片曲线保持实体色彩，形成“描摹字帖”隐喻
•	跨片对齐使用 4.6.5 的情绪进度 u(t)：对当前片与参考片分别计算 u_cur(k)、u_ref(m)，再将参考片曲线 C_ref(m) 重新参数化为 C_ghost(k) = C_ref(u_ref^{-1}(u_cur(k)))；短片/长片自动在“情绪旅程”上对齐
•	UI 主时间轴仍展示线性时间 t 与归一化时间 \tau = t / T_{film}，但 Ghost 曲线的采样点由情绪进度映射决定；无需手调偏移
•	HUD 中同时显示两部片在该情绪进度下的数值、线性时间（分钟）、以及 u(t) 百分比，便于判断谁“情绪旅程领先/滞后”

10.7 技术实现提示与灵感来源
•	优先使用 Swift Charts（macOS 13+）绘制曲线，配置 .interpolationMethod(.catmullRom) 及渐变填充即可满足大部分需求
•	若曲线点数巨大，可降采样或改用 MetalKit Shader 实现 120 fps 绘制
•	参考 UI：Linear（高级渐变）、Activity Monitor（极淡网格）、Blackmagic Cloud/DaVinci（工业按钮）、Oura/Apple Health（把数据当生活方式）
•	所有按钮/面板维持深灰层级，遵循“默认暗色 + 荧光点缀”的 Cinematic Dark Mode