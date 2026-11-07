# GitHub 发布 AMMoire 0.2.0

本目录是独立发布副本，仓库名沿用 `ammoire`，GitHub Pages 关闭。源代码、测试、文献、小型输入和派生基准在发布范围内。本机路径和接力源目录见仓库外的 `DSH_UPLOAD.md`。

## 上传之前

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
python -m ammoire.wannier
python -m ammoire.vasp_bridge --self-test
python -m ammoire demo --shells 2 --points 3 --check-cutoff --out results/release-continuum
python -m ammoire.cairo_demo --points 12 --out results/release-cairo
python -m examples.commensurate_demo --points 8 --out results/release-commensurate
git status --short
git diff --cached --stat
```

本机已验证 Python 3.12.14、NumPy 2.3.5，并完成独立 wheel 构建与源码外安装；具体证据见 `validation/release_checks.json`。远端 Actions 须在推送后查看。

## 已有远端仓库时

先确认用户此前的 `ammoire` URL，并读取其当前分支、历史和用户改动。在新的本地克隆中按发布文件清单合并此副本，保留远端 `.git` 和现有历史。检查差异后创建正常提交，再推送到所选分支。不要把此副本的空历史强制覆盖远端，也不要擅自删除远端独有文件。若有冲突，依据实际变更逐项合并。

## 此前尚未创建远端时

源副本已初始化 `main` 并暂存发布文件，尚未提交和配置远端。使用用户已确认的账户、可见性和 Git 身份创建空仓库；如本机身份未设置，配置该仓库的本地身份后提交。

```sh
git commit -m "Release AMMoire 0.2.0 with commensurate twist controls"
git remote add origin https://github.com/YOUR-ACCOUNT/ammoire.git
git push -u origin main
```

不要将 `YOUR-ACCOUNT` 原样使用。推送后报告仓库 URL、提交 SHA、分支及 CI 状态。

仓库沿用此前未指定再分发许可的状态，许可证由项目所有者另行选择。材料原始输出、授权赝势、波函数、凭据及 dsh 会话资料不在本发布清单中。`benchmark/twist_history` 属于此前计算的派生归档，本轮不把它们标为新的材料重算。
