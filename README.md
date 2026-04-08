# sanger_plot

基于 FastAPI 的 Sanger AB1 文件与参考基因组比对服务，支持单文件、批量分析，以及返回 PNG 图像或 JSON 结果。

## 功能特性

- 单个 AB1 文件比对
- 批量 AB1 文件处理
- 返回 JSON 数据或 PNG 图片
- 使用 BWA 进行序列比对
- 自动检测反向测序
- 支持 Docker / Docker Compose 部署

## 本地运行

### 1. 安装依赖

```bash
pip install -r requirements.txt
```

还需要系统中可用的 `bwa` 命令。

### 2. 配置参考基因组

程序默认读取以下环境变量：

- `DEFAULT_REF_FILE`
- `HG19_REF_FILE`
- `HG38_REF_FILE`
- `OUTPUT_DIR`
- `TEMP_DIR`

示例：

```bash
export HG19_REF_FILE=/data/reference/hg19/hs37d5.fa
export HG38_REF_FILE=/data/reference/hg38/hg38.fa
export DEFAULT_REF_FILE=$HG19_REF_FILE
```

### 3. 启动服务

```bash
python api.py
```

服务默认启动在 [http://localhost:8000](http://localhost:8000)。

## Docker 部署

### 目录约定

建议将参考基因组放在项目下：

```text
data/
  reference/
    hg19/
      hs37d5.fa
      hs37d5.fa.amb
      hs37d5.fa.ann
      hs37d5.fa.bwt
      hs37d5.fa.pac
      hs37d5.fa.sa
    hg38/
      hg38.fa
      hg38.fa.amb
      hg38.fa.ann
      hg38.fa.bwt
      hg38.fa.pac
      hg38.fa.sa
```

如果没有索引文件，应用首次运行时会尝试执行 `bwa index`。

### 使用 Docker Compose

```bash
docker compose up -d --build
```

默认行为：

- 服务端口映射为 `8000:8000`
- 输出图片持久化到 `./outputs`
- 参考基因组以只读方式挂载到容器内 `/data/reference`

停止服务：

```bash
docker compose down
```

### 使用 Docker 命令

构建镜像：

```bash
docker build -t sanger-plot:latest .
```

运行容器：

```bash
docker run -d \
  --name sanger_plot \
  -p 8000:8000 \
  -e DEFAULT_REF_FILE=/data/reference/hg19/hs37d5.fa \
  -e HG19_REF_FILE=/data/reference/hg19/hs37d5.fa \
  -e HG38_REF_FILE=/data/reference/hg38/hg38.fa \
  -v ./outputs:/app/outputs \
  -v ./data/reference:/data/reference:ro \
  sanger-plot:latest
```

## API 文档

启动后可访问：

- Swagger UI: [http://localhost:8000/docs](http://localhost:8000/docs)
- ReDoc: [http://localhost:8000/redoc](http://localhost:8000/redoc)

## API 端点

### `GET /`

返回 API 基本信息和可用参考基因组版本。

### `GET /api/health`

返回服务健康状态与当前检测到的可用参考基因组。

### `GET /full_sequence.html`

返回内置测试页面。

### `POST /api/align`

单个 AB1 文件比对。

表单参数：

- `file`: AB1 文件，必填
- `genome_version`: `hg19` 或 `hg38`，默认 `hg19`
- `window_size`: 显示窗口大小，默认 `10`
- `ref_file`: 自定义参考基因组路径，可选
- `chrom`: 手动指定染色体，可选
- `pos`: 手动指定位置，可选
- `return_json`: 是否返回 JSON，默认 `false`

### `POST /api/align/batch`

批量 AB1 文件比对。

表单参数：

- `files`: 多个 AB1 文件，必填
- `genome_version`: `hg19` 或 `hg38`，默认 `hg19`
- `window_size`: 显示窗口大小，默认 `10`
- `ref_file`: 自定义参考基因组路径，可选
- `return_json`: 是否返回 JSON，默认 `true`
- `chrom`: 手动指定染色体，可选
- `pos`: 手动指定位置，可选

### `POST /api/align/json`

返回 JSON 结果，不生成图片。

## 调用示例

### 单文件返回 JSON

```bash
curl -X POST "http://localhost:8000/api/align" \
  -F "file=@sample.ab1" \
  -F "genome_version=hg19" \
  -F "window_size=10" \
  -F "return_json=true"
```

### 批量处理

```bash
curl -X POST "http://localhost:8000/api/align/batch" \
  -F "files=@sample1.ab1" \
  -F "files=@sample2.ab1" \
  -F "genome_version=hg19" \
  -F "return_json=true"
```

## 输出说明

- 生成的 PNG 文件会保存到 `OUTPUT_DIR`
- API 同时会暴露 `/outputs/<filename>` 静态路径访问生成图片
- 批量接口返回的 `image_url` 会自动使用当前请求地址生成

## 依赖项

- FastAPI
- Uvicorn
- Biopython
- NumPy
- Matplotlib
- BWA

## 许可证

MIT License
