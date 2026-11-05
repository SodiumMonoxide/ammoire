# GitHub 发布

从本目录发布。它是 AMMoire 的独立仓库，包含源码、测试、文献、可分享的 VASP 输入和派生基准。原始 VASP 输出、本机计算中间文件、授权赝势、旧日报库及 dsh 会话资料留在本地研究目录。

创建一个空的 GitHub 仓库，选择可见性，然后在本目录核对待提交文件。Git 用户身份按需配置为该仓库的本地配置。仓库目前没有指定再分发许可，由项目所有者选择是否及如何授权再使用。

```sh
git status --short
git diff --cached --stat
git config user.name "YOUR NAME"
git config user.email "YOUR EMAIL"
git add -A
git commit -m "Initial AMMoire release"
git remote add origin https://github.com/YOUR-ACCOUNT/YOUR-REPO.git
git push -u origin main
```

推送后可在 GitHub 的 Actions 查看 CI。随后在本地 dsh 客户端打开克隆目录，使用 [HANDOFF.md](HANDOFF.md) 顶部的接力任务。本机复核已有 VASP 原始结果时，将其所在目录作为 `analyse.py` 的输入，不必复制到仓库。
