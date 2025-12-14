# ArcScope

> 设计文档（建议先读）：`Arcscope.docx` / `arcscope.md` / `ArcScope_CoreContract.md`

## 一句话定义

**ArcScope 是一个完全本地运行的电影结构分析系统**：

它把一部影片统一转换为 **1Hz 的多模态时间序列**，通过**稳健统计与 PCA** 自动建模节奏、声音、色彩、信息与表演，并在同一时间轴上生成**结构曲线、问题诊断与多尺度叙事分段**，用于复盘、学习与创作校准。

这句话隐含了所有核心设计原则：
- ✅ 本地运行（无云端依赖）
- ✅ 统一时间轴（1Hz）
- ✅ 无参数（自适应）
- ✅ PCA 自动权重
- ✅ 曲线 + 诊断 + 结构

---

## 项目概述（来自 1.0 设计文档）

ArcScope 面向电影创作者（导演、剪辑、制片、研究者）的本地电影结构分析工具：把一部影片转换为若干条随时间变化的“结构曲线”，帮助你看清节奏、情绪、信息与视听关系。

主要使用场景：
- 剪辑/混音阶段的复盘与调节
- 对经典影片的结构拆解学习
- 自己作品的二刷分析与问题段定位

关键特性：
- 完全本地分析（无云端上传，规避版权与隐私问题）
- 统一时间轴：所有曲线以秒为横坐标，归一化到 `[0,1]` 便于叠加对比
- 多模态综合：剪辑节奏、运动、声音、色彩、字幕信息、表演情绪
- 自动诊断：标出“拖沓 / 情绪空转 / 视听错位 / 解释压死情绪”等区段
- 参考对比：叠加参考影片曲线做结构对照

---

## 功能清单（1.0 与预留）

### 1.0 必须具备（核心底座）

- 影片导入与分析启动（本地文件）
- 曲线（统一时间轴）：
  - `Pace(t)`：节奏（剪辑 + 运动）
  - `Sound(t)`：声音能量/压力
  - `Color(t)`：色彩情绪/视觉温度
  - `Info(t)`：信息密度（语言/视觉/事件负荷）
  - `Arousal(t)`：综合唤醒
  - 可选：`FaceAffect(t)`：表演存在度与情绪强度（可用则启用）
- 自动诊断区段（IssueSegments）：至少覆盖低活跃、情绪空转、视听错位、信息压制等类型
- 结构分段（StructureSegments）：`Shot → Scene → Sequence → Act`
- 本地影片库（SQLite 持久化与浏览）
- 参考库：将影片标记为参考，并在视图中叠加 ghost 曲线
- 导出：海报（PNG/PDF）与标准报告（PDF/Markdown）

### 预留增强（不阻塞 1.0，但要留接口）

- 以 CoreML 多模态模型替代部分 arousal 合成
- 文本情绪细分（字幕 valence/arousal）
- 风格/导演地图（PCA/UMAP 空间可视化）

---

## 系统架构总览

ArcScope 采用**严格三层架构**，实现平台无关的核心算法与 macOS 平台的无缝集成：

```
┌─────────────────────────────────────────┐
│  App Layer (SwiftUI)                    │
│  • UI渲染与交互                          │
│  • 曲线可视化、热力图、问题列表           │
└──────────────┬──────────────────────────┘
               │ C接口调用
┌──────────────▼──────────────────────────┐
│  Bridge Layer (Objective-C++)           │
│  • 视频解码与特征提取                     │
│  • Apple框架集成（Vision、AVFoundation）  │
└──────────────┬──────────────────────────┘
               │ FeatureSample[]
┌──────────────▼──────────────────────────┐
│  Core Layer (C++17)                     │
│  • 纯算法引擎：特征融合→曲线生成→诊断      │
│  • FilmEngine编排整个分析流程             │
└──────────────┬──────────────────────────┘
               │
               ▼ SQLite
          [数据库持久化]
```

### 总体原则（非常重要）

| 层级 | 职责 |
|------|------|
| **Bridge** | 只"观测"，不"判断" |
| **Core** | 只"建模"，不"采集" |
| **App** | 只"呈现"，不"计算" |
| **SQLite** | 是唯一真相源 |

---

## 端到端数据流（冻结版）

```
video file
  ↓
[ Bridge Layer ]
  ↓  ArcScopeBridgeSample[0..N-1]   (1Hz, raw observations)
[ Core Layer ]
  ↓  FeatureSeries / ZSeries
  ↓  Curves (0~1)
  ↓  IssueSegments / StructureSegments
[ SQLite ]
  ↓
[ SwiftUI App ]
```

**任何模块不允许跨越这一方向。**

### 详细数据流

#### 阶段 1：原始数据提取 (Bridge Layer)

```
video.mp4 输入
    │
    ├→ VideoFeatureExtractor (AVFoundation)
    │   • 720p@10fps 解码
    │   • 光流计算 → motion_amplitude
    │   • 镜头检测 → cut_density
    │
    ├→ AudioAnalysisEngine (CoreAudio)
    │   • PCM 提取 → audio_rms, transients
    │   • 频谱分析 → spectral_balance
    │
    ├→ FaceTrackingEngine (Vision + CoreML)
    │   • 人脸检测 + 多目标跟踪
    │   • 表情识别 → face_valence, face_arousal
    │
    ├→ ColorSciencePipeline
    │   • RGB → ACEScg → CAM16 色彩空间
    │   • brightness, saturation, warmth, grayness
    │
    └→ NLPProcessor (NaturalLanguage)
        • 字幕解析 → 分词、词性标注
        • word_count, concept_density

输出：ArcScopeBridgeSample[N]
      每秒 18 个原始观测值
```

#### 阶段 2：特征扩展 (Core Layer)

```
FeatureSample[] 输入
    │
    ├→ SoundAnalyzer
    │   • Loudness, RhythmEnergy, Chroma
    │   • Key detection, Harmonic tension
    │
    ├→ ColorAnalyzer
    │   • Color harmony, Warmth evolution
    │   • Saturation dynamics
    │
    ├→ FaceAnalyzer
    │   • 主角识别（自适应算法）
    │   • FacePresence × ExpressionIntensity
    │
    ├→ MotionAnalyzer
    │   • ASL (平均镜头长度)
    │   • Shot segmentation
    │
    └→ InfoAnalyzer
        • VerbalLoad, VisualLoad
        • Exposition probability

输出：AudioFeatures, ColorFeatures, etc.
```

#### 阶段 3：曲线生成 (CurveEngine)

使用 **PCA** 自动权重分配，生成 6 条核心曲线：

```
CurveEngine (PCA-based synthesis)
    │
    ├→ Pace = PCA1([-ASL, CutDensity, Motion])
    │   「节奏：剪辑密度 + 运动强度」
    │
    ├→ Sound = PCA1([Loudness, Rhythm, Harmony])
    │   「音频：响度 + 节奏 + 和声」
    │
    ├→ Color = PCA1([Warmth, Brightness, Saturation])
    │   「色彩：暖度 + 亮度 + 饱和度」
    │
    ├→ Info = PCA1([VerbalLoad, VisualLoad, EventLoad])
    │   「信息：对话 + 文本 + 事件」
    │
    ├→ Face = FacePresence × ExpressionIntensity
    │   「表演：人脸存在度 × 情绪强度」
    │
    └→ Arousal = PCA1([Pace, Sound, Color, Info, Face])
        「唤醒度：所有维度的综合」

所有曲线经过：
  1. Robust Z-score (中位数 + MAD)
  2. Sigmoid normalization → (0, 1)
```

#### 阶段 4：自动诊断 (DiagnosticsEngine)

基于**自适应阈值**（分位数）检测 13 种问题类型：

| 问题类型 | 检测逻辑 | 示例 |
|---------|---------|------|
| **LowActivity** | Pace < Q_P(0.25) 持续过长 | 节奏低迷片段 |
| **OvercutFlat** | Pace > Q_P(0.75) BUT Arousal < Q_A(0.25) | 剪辑频繁但无情绪 |
| **AudioVisualMisalign** | Sound 和 (Pace\|Color) 反相关 | 画面安静但音乐激烈 |
| **InfoKillsEmotion** | Info high + Arousal drops | 对白过多压制情绪 |
| **HighTensionLowPace** | Sound/Color high + Pace low | 视觉张力但节奏慢 |

严重度计算：`∫(偏离量) dt` （面积积分，无固定阈值）

#### 阶段 5：结构分割 (StructureAnalyzer)

```
层次化分割（4个层级）：

Shot (镜头)
  ↓ 光流 + 硬切检测
Scene (场景)
  ↓ Changepoint detection on [Pace, Sound, Info]
  ↓ Scene typing: LowEnergy, BuildUp, Peak, etc.
Sequence (序列)
  ↓ Greedy agglomeration (目标 4-12 个)
Act (幕)
  ↓ Hierarchical split (2-4 幕)
```

每个 Scene 附带：
- **类型标签**（自动推断）
- **色彩指纹**（主导色 + 饱和度分布）

---

## Bridge Layer（Objective-C++）最终职责与契约

### Bridge 的唯一职责

**把影片转换成"每秒一行的原始观测表"。**

### Bridge 允许做的

- ✅ 解码（AVFoundation）
- ✅ 下采样（720p / 10fps）
- ✅ 聚合（10fps → 1Hz）
- ✅ 调用 Apple 专有能力（Vision / CoreML / Accelerate / NL）

### Bridge 严禁做的

- ❌ PCA
- ❌ 阈值判断（"高/低/异常"）
- ❌ 曲线合成
- ❌ 诊断 / 结构判断
- ❌ "主角是谁"的最终裁决（只提供原始轨迹）

### ArcScopeBridgeSample（冻结字段定义）

这是 Core 与 Bridge 的硬契约，1.0 不再随意加字段。

**每秒一条，长度 = ⌈runtime⌉**

| 分类 | 字段 | 说明 |
|------|------|------|
| **时间** | `time_seconds` | k + 0.5 |
| **剪辑/运动** | `cut_density` | 该秒内 cut 数 |
| | `motion_amplitude` | 光流幅值均值 |
| **音频** | `audio_rms` | RMS |
| | `audio_transient` | 瞬态强度 |
| | `spectral_balance` | 高频/低频比 |
| **色彩** | `brightness` | 感知亮度 |
| | `saturation` | 饱和度 |
| | `warmth` | 暖冷 |
| | `grayness` | 去色程度 |
| | `avg_hue` | 主色调 |
| **文本** | `subtitle_density` | 字/秒 |
| | `exposition_probability` | 解释性 |
| **信息** | `info_density` | Bridge 级粗估 |
| **人脸** | `face_presence` | 任一脸是否存在 |
| | `face_arousal` | 模型输出 |
| | `face_valence` | 模型输出 |
| | `face_confidence` | 检测置信度 |

**注意**：
- ASL / Chroma / Key / VerbalLoad / VisualLoad **永远不在 BridgeSample 中**
- BridgeSample 是"测量"，不是"理解"

### C 接口（供 Swift 调用）

```c
// 端到端分析（视频 → SQLite）
arcscope_analyze_film(video_path, db_path, film_id, title, year, director,
                      progress_cb, user_data)

// 加载分析结果
arcscope_load_film_detail(db_path, film_id) → ArcScopeFilmDetail*

// 列举影片库
arcscope_list_films(db_path) → ArcScopeFilmList*

// 直接驱动分析（用于测试）
arcscope_bridge_analyze(samples[], n_samples, ...)
```

---

## Core Layer（C++17）最终架构（冻结）

### Core 的一句话职责

**把 1Hz 的原始观测转化为"自适应、可解释的结构曲线与区段"。**

### Core 的模块边界（最终）

```
FilmEngine
 ├─ FeatureExpansion
 │   ├─ SoundAnalyzer
 │   ├─ MotionAnalyzer
 │   ├─ ColorAnalyzer
 │   ├─ InfoAnalyzer
 │   └─ FaceAnalyzer
 ├─ Normalization (robust)
 ├─ CurveEngine
 ├─ ArousalEngine
 ├─ DiagnosticsEngine
 ├─ StructureAnalyzer
 └─ SQLiteStore
```

### 1️⃣ FeatureExpansion（从 raw 到语义特征）

**输入**：BridgeSample[]
**输出**：FeatureSeries（仍是 1Hz）

| Analyzer | 输出特征 |
|----------|---------|
| **SoundAnalyzer** | Loudness / RhythmEnergy / Chroma / Key / HarmonicTension |
| **MotionAnalyzer** | ASL / CutDensity / Motion |
| **ColorAnalyzer** | Warmth / Brightness / Saturation / Harmony |
| **InfoAnalyzer** | VerbalLoad / VisualLoad / EventLoad |
| **FaceAnalyzer** | FaceAffect[k] / FacePresenceMask[k] |

**FaceAnalyzer（关键）**：
- 全片统计 → dominant face 决策
- 若失败：`face_enabled = false`
- **只有 Core 有资格决定"是否存在主角"**

### 2️⃣ Normalization（全系统统一）

**唯一合法公式**：

```
Z(k) = (X_k - median(X)) / (MAD(X) + ε)
X̃(k) = σ(Z(k))  // sigmoid
```

- 无分支
- 无 if
- 无手调参数

### 3️⃣ CurveEngine（PCA 主轴建模）

| 曲线 | 输入（Z） |
|------|---------|
| **Pace** | [-ASL, CutDensity, Motion] |
| **Sound** | [Loudness, Rhythm, Tension] |
| **Color** | [Warmth, Brightness, Saturation] |
| **Info** | [Verbal, Visual, Event] |
| **Face** | FacePresence × Expression |
| **Arousal** | [Pace, Sound, Color, Info, Face?] |

**规则冻结**：
1. PCA PC1
2. 若与直觉方向相反 → 乘 -1
3. Face 不启用 → 从 PCA 维度中移除

### 4️⃣ DiagnosticsEngine（问题检测）

**统一逻辑**：
- 阈值 = 分位数
- 持续时间 = 稳定段中位长度
- 严重度 = 偏离面积积分

`IssueSegment` 是 Core 的最终语义输出之一。

### 5️⃣ StructureAnalyzer（结构分段）

**层级冻结**：

```
Shot  →  Scene  →  Sequence  →  Act
```

- **Scene**：Changepoint on [Pace, Sound, Info]
- **Sequence**：相邻 Scene 聚合
- **Act**：宏观层级（2–4）

**存储规则冻结**：
- 所有结构段写入 `segments`
- `type` 使用前缀：
  - `issue:low_activity`
  - `scene:build_up`
  - `sequence:conflict`
  - `act:II`

### 核心模块文件

| 模块 | 文件 | 职责 |
|------|------|------|
| **FilmEngine** | `core/include/arcscope/FilmEngine.h` | 主编排器：串联 5 个分析阶段 |
| **CurveEngine** | `core/include/arcscope/CurveEngine.h` | PCA-based 曲线生成 |
| **DiagnosticsEngine** | `core/include/arcscope/Diagnostics.h` | 自适应诊断（13 种 Issue） |
| **StructureAnalyzer** | 内嵌于 FilmEngine | 多尺度结构分割 |
| **SQLiteStore** | `core/include/arcscope/SQLiteStore.h` | 数据库持久化（仅存储） |
| **Normalization** | `core/include/arcscope/Normalization.h` | 稳健统计工具库 |

**关键特性**：
- ✅ **零 Apple 依赖**（可移植到 Windows/Linux）
- ✅ **无硬编码阈值**（所有参数来自数据分位数）
- ✅ **PCA 自动权重**（无需手调）

---

## SQLite（最终口径）

**SQLite 是分析结果的唯一事实源。**

### 数据库模式

```sql
-- Core tables
films(id INTEGER PK, title, year, director, runtime_sec, ...)
curves(id, film_id FK, curve_type TEXT, fps=1.0, length, data_blob BLOB)
metrics(id, film_id FK, name TEXT, value REAL)
segments(id, film_id FK, type TEXT, start_sec, end_sec, severity)

-- Invariants
-- R3: Every curves.fps = 1.0
-- R8: Every curves.data_blob = float32 array (binary)
-- Time: times[i] = i + 0.5 for i = 0..length-1
```

### 存储规则

| 数据类型 | 规则 |
|---------|------|
| **Curves** | 只存 0~1 曲线 |
| **Face** | • face_presence 永远可存<br>• face_affect 只有启用时才存 |
| **Metrics** | • PCA loadings<br>• 分位数阈值 |
| **Segments** | • Issue + Structure 共用表<br>• 用 type 前缀区分 |

---

## SwiftUI App（最佳实践）

### App 的原则

| 原则 | 说明 |
|------|------|
| **零算法** | 不做 PCA、不做特征提取 |
| **零判断** | 不判断阈值、不检测问题 |
| **零 PCA** | 不计算主成分 |
| **只渲染** | 从 SQLite 拉数据 → 渲染 |
| **只交互** | 游标、hover、时间跳转 |

### ViewModel 唯一职责

1. 从 SQLite 拉数据
2. 做时间 → x 映射
3. reference 重采样（情绪进度 u(t)，可延后）

### UI 组件

| 组件 | 文件 | 功能 |
|------|------|------|
| **ArcScopeViewModel** | `app/Sources/.../ViewModels/ArcScopeViewModel.swift` | 应用状态管理 |
| **BridgeBackend** | `app/Sources/.../Backend/BridgeBackend.swift` | Bridge 接口包装 |
| **CurveTimelineView** | `app/Sources/.../Views/CurveTimelineView.swift` | 6 条曲线渲染 + 交互 |
| **DiagnosticsHeatmap** | `app/Sources/.../Views/DiagnosticsHeatmap.swift` | 问题热力图 |
| **SceneBandsView** | `app/Sources/.../Views/SceneBandsView.swift` | 4 层结构色带 |
| **IssueListView** | `app/Sources/.../Views/IssueListView.swift` | 问题列表 + 时间跳转 |

**UI 亮点**：
- 🎨 6 条曲线叠加渲染（面积图 + 辉光效果）
- 🖱️ 交互游标 + HUD（悬停显示数值）
- 🔥 热力图（颜色=类型，深度=严重度）
- 📊 多层色带（Shot/Scene/Sequence/Act）

---

## ArcScope 1.0 冻结清单（非常重要）

### ✅ 1.0 必须完成

- [x] 三层严格分离
- [x] 1Hz 时间轴
- [x] 6 条曲线（Pace / Sound / Color / Info / Face / Arousal）
- [ ] Face 可禁用逻辑
- [ ] 至少 5 类 Issue
- [ ] Scene / Sequence / Act 分段
- [ ] SQLite 完整持久化

### ❌ 明确延后（写进 TODO）

- [ ] 多模态深度模型
- [ ] UMAP / 风格地图
- [ ] 自动报告 NLP 总结
- [ ] 云端 / 协作

---

## 构建与运行

### 前置要求

- macOS 13+ (Apple Silicon 或 Intel)
- Xcode 14+ (Command Line Tools)
- CMake 3.20+
- Swift 5.7+

### 构建 FilmEngine (Core Layer)

```bash
CXX=/usr/bin/clang++ cmake -S core -B core/build
cmake --build core/build
./core/build/arcscope_demo   # 生成 arcscope.db（合成数据）
```

CLI demo 运行完整流程：normalization → PCA curves → diagnostics → SQLite persistence。

### 构建 Bridge Layer

```bash
CXX=/usr/bin/clang++ cmake -S bridge -B bridge/build
cmake --build bridge/build
./bridge/build/arcscope_bridge_demo   # 调用 bridge API 端到端
```

Bridge 暴露 `arcscope_bridge.h`，Swift/Objective-C 可通过系统模块 `CAArcScopeBridge` 导入。

### 真实视频分析

```bash
./bridge/build/arcscope_video_demo /path/to/movie.mp4
```

使用 AVFoundation 解码视频，下采样至 720p@10fps，计算光流、人脸、色彩、音频特征，生成 1Hz 样本流，调用 FilmEngine 分析，结果写入 `video_demo.db`。

### 构建 SwiftUI App

```bash
# 先构建核心与桥接库（仅需一次）
CXX=/usr/bin/clang++ cmake -S core -B core/build && cmake --build core/build
CXX=/usr/bin/clang++ cmake -S bridge -B bridge/build && cmake --build bridge/build

# 分析任意影片，写入 SQLite
./bridge/build/arcscope_video_demo /path/to/movie.mp4

# 启动 SwiftUI 应用（会链接 libarcscope_bridge.a）
cd app
swift build  # requires macOS 13+ SDK
swift run ArcScopeApp
```

运行与数据位置说明：
- 若未选择数据库，App 默认把数据库写到 `~/Library/Application Support/ArcScope/ArcScope.sqlite`
- 通过 UI 的文件选择器选择数据库后，会保存 **security-scoped bookmark**，重启后仍可访问（为 App Sandbox / 上架准备）
- 视频文件同样会保存 bookmark（用于导出海报的关键帧缩略图读取）

> **Sandbox note**：如果你在沙盒环境（App Sandbox）运行，所有文件访问都应通过系统文件选择器/书签访问模型完成；不要依赖固定路径（如 `~/Desktop/...`）。

---

## macOS 菜单与快捷键（HIG 基建）

- 菜单项：选择数据库、分析新影片、参考库、导出、设置（标准 Settings 窗口）
- 默认快捷键：`⌘O` 选择数据库；`⌘N` 分析新影片；`⌘R` 参考库；`⌘E` 导出（需要已有分析结果）

---

## 本地化（Localization）

本项目已将 UI/错误文案等收口到本地化资源中：

- 资源：`app/Sources/ArcScopeApp/Resources/*.lproj/Localizable.strings`
- SwiftPM 默认语言：`app/Package.swift` 中 `defaultLocalization`

新增/修改文案的约束：
- 所有用户可见字符串必须走 key（避免硬编码）
- 时间/百分比/小数等展示必须使用系统格式化（避免地区格式翻车）

---

## 数据与隐私

- ArcScope 的分析与渲染完全在本机执行，不依赖云端服务。
- 分析结果持久化为 SQLite（本地文件），可由用户自行删除数据库文件来清理数据。

---

## App Store / Sandbox 准备要点（简版）

如果你要做可上架的 macOS app（App Sandbox 开启）：
- 文件访问：始终通过 `fileImporter`/`NSOpenPanel` + security-scoped bookmark
- 默认存储：使用 Application Support（不要写 Desktop/Downloads）
- 权限：Entitlements 最小化

完整验收清单见：`docs/AppleGuidelinesAcceptanceChecklist.md`

### 集成到 Xcode

1. 将 `bridge/build/libarcscope_bridge.a` 添加到 Xcode 项目
2. 通过 bridging header 暴露 `arcscope_bridge.h`
3. 实现 `ArcScopeBackend` 包装 `arcscope_bridge_*` 函数
4. 将 `ArcScopeFilmDetail` 转换为 Swift `FilmAnalysis` 模型
5. 用 bridge-backed 实现替换 `MockBackend`

---

## 仓库结构

```
Arcscope/
├── core/                    # FilmEngine C++17 核心引擎
│   ├── include/arcscope/
│   │   ├── FilmEngine.h           # 主编排器
│   │   ├── FilmTypes.h            # 数据结构
│   │   ├── CurveEngine.h          # PCA 曲线生成
│   │   ├── Diagnostics.h          # 问题检测
│   │   ├── SQLiteStore.h          # 持久化
│   │   ├── Normalization.h        # 稳健统计
│   │   └── features/              # 特征模块
│   └── src/                       # 实现文件
│
├── bridge/                  # Objective-C++ 桥接层
│   ├── include/
│   │   └── arcscope_bridge.h      # C 接口
│   ├── src/
│   │   └── arcscope_bridge.mm     # Bridge 实现
│   ├── VideoFeatureExtractor.h/mm # 视频解码
│   ├── AudioAnalysisEngine.h/mm   # 音频分析
│   ├── FaceTrackingEngine.h/mm    # 人脸检测
│   ├── ColorSciencePipeline.h/mm  # 色彩转换
│   └── tests/                     # 测试程序
│
├── app/                     # SwiftUI macOS 应用
│   ├── Package.swift
│   └── Sources/ArcScopeApp/
│       ├── Backend/               # Bridge 包装
│       ├── Models/                # 数据模型
│       ├── ViewModels/            # 状态管理
│       └── Views/                 # UI 组件
│
├── db/                      # 数据库目录
├── arcscope.md              # 完整设计文档
├── ArcScope_CoreContract.md # 架构约束
└── README.md                # 本文件
```

---

## 核心技术栈

| 层级 | 语言 | 关键库/框架 | 代码量 |
|------|------|-----------|--------|
| **Core** | C++17 | SQLite3, 标准库 | ~1724 行 |
| **Bridge** | Objective-C++ | AVFoundation, Vision, CoreML, Accelerate, NaturalLanguage | ~800 行 |
| **App** | Swift | SwiftUI, AppKit | ~600 行 |

---

## 架构约束（不可违反）

详见 `ArcScope_CoreContract.md`，以下为关键约束：

| 约束 | 要求 |
|------|------|
| **R1: 1Hz时间轴** | Core仅操作1Hz样本 → sample k代表[k,k+1) |
| **R2: Bridge输出** | 仅原始观测，不融合/规范化/分类 |
| **R3: 数据库fps** | 所有curves.fps = 1.0（强制不变量） |
| **R4: 语义边界** | Bridge≠Core≠Store（严格分离） |
| **R5: 执法规则** | 违反→代码评审时BLOCK合并 |
| **R6: ID二重性** | film_logical_id（用户级）+ film_pk（数据库级） |
| **R7: 内存所有权** | 谁分配谁释放（bridge_free_*） |
| **R8: Blob格式** | curves.data_blob = float32[]（二进制） |

**违反任何约束都会在 code review 时被拒绝。**

---

## 扩展路线图

1. **媒体提取完善**：集成 AVFoundation + Vision 到 Bridge，实现完整的 10fps → 1Hz 特征流
2. **诊断丰富**：扩展 DiagnosticsEngine 到 13 种 Issue 类型
3. **SQLite 同步**：实现 Swift 端查询/读取模型，支持影片库浏览
4. **参考对比**：添加 u(t) 驱动的 ghost 曲线 + PCA 对齐
5. **导出功能**：支持海报/PDF 报告导出

---

## License

本项目遵循设计文档 `arcscope.md` 中的规范，用于电影分析与创作辅助。

---

**最后更新**: 2025-12-13
