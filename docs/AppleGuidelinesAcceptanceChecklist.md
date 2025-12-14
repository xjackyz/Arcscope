# ArcScope：Apple 官方“总纲”对照检查与 1.0 验收清单

> 目标：让 ArcScope 在 **macOS 交互一致性 / 本地化（含混合语言内容）/ 沙盒与文件访问 / 隐私表达 / 上架合规** 上，能按 Apple 的官方总纲逐项验收、可测、可复盘。

## 0. 官方“总纲”（你需要逐条过一遍的 Apple 文档）

这些是“总纲级”文档：它们不是教你某个 API，而是规定 **什么算像 Mac、什么算合规、什么算本地化做对**。

- Human Interface Guidelines（HIG，macOS）：窗口与层级、布局密度、菜单栏与快捷键、可预测的窗口行为、避免把关键控件放在窗口底部等。
- Apple Localization 总页：本地化流程、格式化（日期/数字/货币/单位）、以及“用户生成内容可以混合多语言”的原则。
- Xcode Localization / String Catalog（`.xcstrings`）：Xcode 15+ 推荐的字符串资产（复数/变体/翻译注释）。
- SwiftUI：Preparing views for localization：SwiftUI 写法如何让可导出字符串、带 hint/comment、可适配文本扩展。
- App Sandbox（macOS）：本地文件型 app 最易踩坑的文件访问/权限模型（用户选取、书签、安全作用域 URL）。
- App Review Guidelines（macOS 相关）：沙盒、文件系统规范、打包提交流程等（尤其是“不要绕过系统文件访问模型”）。

## 1) 多语言（ArcScope 这种“重 UI + 重术语”）落地方式

### 1.1 文本体系：统一入口，避免“手搓本地化”

- [ ] **所有用户可见字符串** 统一收口（按钮、提示、错误、空状态、设置项、导出报告标题/字段名、Issue 名称与解释等）。
- [ ] **同一个英文词在不同语境** 必须有不同 key 或提供明确 hint/comment（避免翻译歧义）。
- [ ] **复数与数量句子**（如“找到 3 段问题”）使用 String Catalog 的复数/变体能力（不要用 if/拼接）。
- [ ] **术语表**：建立术语一致性（比如 “Shot/Scene/Sequence/Act/Issue” 以及曲线名、诊断名）。

> 备注：当前仓库是 Swift Package（`app/Package.swift`），是否能“完整享受”`.xcstrings` 取决于你最终的构建与打包方式（Xcode 工程/Archive）。如果最终走 Xcode 上架，建议 UI 文本最终迁移到 `.xcstrings`；否则至少要用 `Resources/*.lproj/Localizable.strings` 做资源本地化。

### 1.2 格式化：禁止字符串拼接

- [ ] 时间/时长：用 `DateComponentsFormatter`（`.positional`，支持 `mm:ss` / `hh:mm:ss`），不要 `String(format: "%02d:%02d")`。
- [ ] 百分比：用 `NumberFormatter`（`.percent`），不要 `String(format: "%.0f%%")`。
- [ ] 数字/小数：用 `NumberFormatter`/`FormatStyle`（避免小数点、分组符号随地区变化导致显示错误）。
- [ ] 单位（分钟、秒、大小等）：用 `MeasurementFormatter` 或系统 `FormatStyle`，不要拼 `m/s`。

### 1.3 文本扩展与布局（多语言压力测试必做）

- [ ] 语言压力测试：`zh-Hans`、`en`、`de`（长文本）至少三组。
- [ ] 地区压力测试：例如 `en_US` vs `fr_FR`（数字分组、小数点、百分号位置等）。
- [ ] 关键 UI 不依赖固定宽度（按钮/标题区能扩展、允许换行策略明确）。
- [ ] 避免把关键按钮放在窗口最底部（窗口下边缘可能被拖出屏幕）。

### 1.4 “内容语言”与“UI 语言”分离（ArcScope 的关键）

- [ ] **数据库内容**（片名、导演、备注、标签、字幕统计等）按原始内容存储与展示，不随 UI 语言强制转写/替换。
- [ ] **Issue 的解释文本** 可以本地化，但 **影片原始元数据** 不应被“本地化覆盖”。
- [ ] 支持 **同一屏混合语言内容**（例如 UI 中文，影片英文标题+日文字幕片段）。

## 2) “像 Mac”且可测的 macOS 体验检查（HIG 落地）

### 2.1 菜单栏与快捷键（macOS 的基本盘）

- [ ] 菜单栏提供关键操作：导入/打开、导出、设置、帮助。
- [ ] 快捷键覆盖高频操作：打开数据库/导入视频、搜索、导出、切换视图/曲线、显示设置等。
- [ ] 设置入口遵循系统习惯（系统 Settings 菜单项/独立 Settings 窗口，而不是“只在某个页面角落放齿轮”）。

### 2.2 窗口与状态可预测

- [ ] 分析进行时，窗口状态/进度展示明确且不遮挡关键控件。
- [ ] 多窗口/单窗口策略清晰（是否允许打开多个影片？导出/设置是否独立窗口？）
- [ ] 默认窗口尺寸与最小尺寸合理（避免小屏幕被迫横向滚动）。

### 2.3 可访问性与系统一致性（别做“只能看不能用”的专业 UI）

- [ ] 不强制全局 `.foregroundStyle(.white)`（需要支持浅色/深色、对比度、系统颜色语义）。
- [ ] 图标/颜色有非颜色语义提示（Issue 不能只靠颜色区分）。
- [ ] 文字可被 VoiceOver 读懂（图标、曲线名、Issue 名称等需要可读 label/hint）。

## 3) 本地文件访问与权限（ArcScope 的高危区）

### 3.1 沙盒最小权限

- [ ] Entitlements 只申请必要能力（先按“最小可用”配置）。
- [ ] 不依赖固定路径（如 `~/Desktop/...`）作为默认数据源或默认数据库位置。

### 3.2 文件访问走系统模型（用户选取 + 书签）

- [ ] 打开视频/数据库使用系统选择器（SwiftUI `.fileImporter` / AppKit `NSOpenPanel`）。
- [ ] **需要跨重启持续访问** 的路径（数据库、长期素材库、导出目录）使用 **security-scoped bookmark** 持久化。
- [ ] 访问 security-scoped URL 时，严格 `startAccessingSecurityScopedResource()` / `stopAccessing...` 配对。
- [ ] 诊断 sandbox 违规：跑一遍 Apple 提供的沙盒违规诊断路径，把告警清零。

## 4) 隐私与信任感（全本地是优势，但要“讲清楚”）

- [ ] 在合适的位置明确：是否联网、是否上传、数据保存在哪里、如何删除数据。
- [ ] 隐私表达不藏在角落一行小字，且与实际行为一致（例如导出、分析是否触发网络）。

## 5) 基于当前仓库的快速现状扫描（可直接开工的 TODO）

### 5.1 本地化与格式化“必改点”

- [ ] 硬编码 UI 字符串大量存在（示例）：
  - `app/Sources/ArcScopeApp/Views/LibraryView.swift:36`
  - `app/Sources/ArcScopeApp/Views/SettingsView.swift:11`
  - `app/Sources/ArcScopeApp/Views/Export/ReportExportView.swift:33`
- [ ] 错误信息当前是单一中文文案（需要进入本地化体系）：
  - `app/Sources/ArcScopeApp/Backend/ArcScopeBackend.swift:15`
  - `app/Sources/ArcScopeApp/Backend/VideoAnalysisService.swift:11`
- [ ] 时间/百分比/单位存在手工拼接与 `String(format:)`（地区格式会翻车）：
  - `app/Sources/ArcScopeApp/Models/TimelineFormatting.swift:4`
  - `app/Sources/ArcScopeApp/Views/LibraryView.swift:188`
  - `app/Sources/ArcScopeApp/Views/Export/ReportExportView.swift:231`
  - `app/Sources/ArcScopeApp/Views/StructureBandsView.swift:136`

### 5.2 沙盒与文件访问“高风险点”

- [ ] 默认数据库硬编码到桌面路径（沙盒 app 默认不可访问 Desktop；且不可复现）：
  - `app/Sources/ArcScopeApp/ViewModels/ArcScopeViewModel.swift:24`
  - `app/Sources/ArcScopeApp/Backend/BridgeBackend.swift:25`
- [ ] `.fileImporter` 已在用（很好），但未见对“跨重启访问数据库”的 bookmark 持久化逻辑（上架前必须补齐）：
  - `app/Sources/ArcScopeApp/Views/LibraryView.swift:65`

### 5.3 HIG（macOS）一致性“缺口”

- [ ] `ArcScopeApp` 目前未定义菜单栏命令与快捷键（建议补 `Commands` / `Settings` scene）：
  - `app/Sources/ArcScopeApp/ArcScopeApp.swift:7`
- [ ] 全局强制白色前景（浅色模式/对比度/可访问性风险）：
  - `app/Sources/ArcScopeApp/Views/LibraryView.swift:30`

## 6) 1.0 验收建议顺序（按风险与收益排序）

1. [ ] **沙盒文件访问模型**：去掉默认 `~/Desktop/...`，改成“首次启动引导选择数据库/目录 + bookmark 持久化 + 访问配对”。
2. [ ] **文本体系**：先把错误/空状态/按钮/设置项统一收口（准备迁移 `.xcstrings` 或 `.strings`）。
3. [ ] **格式化改造**：时间/百分比/单位全部换系统 formatter（同时做 `de` 压力测试）。
4. [ ] **HIG 基建**：菜单栏/快捷键/Settings 窗口/帮助入口齐全。
5. [ ] **隐私表达**：写清楚“全本地、不上传”的范围与例外（如未来可能的在线功能）。

