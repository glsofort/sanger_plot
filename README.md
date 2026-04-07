# sanger_plot

# 基因组比对API

基于 FastAPI 的 AB1 文件与参考基因组比对服务。

## 功能特性

- 单个AB1文件比对
- 批量AB1文件处理
- 返回JSON数据或PNG图片
- 使用BWA进行序列比对
- 自动检测反向测序

## 安装依赖

```bash
pip install -r requirements.txt
```

## 启动服务

```bash
python api.py
```

服务将在 `http://localhost:8000` 启动。

## API文档

启动服务后，访问以下地址查看自动生成的API文档：

- Swagger UI: `http://localhost:8000/docs`
- ReDoc: `http://localhost:8000/redoc`

## API端点

### 1. 根路径
```
GET /
```
返回API基本信息。

### 2. 健康检查
```
GET /api/health
```
返回服务健康状态。

### 3. 单个文件比对
```
POST /api/align
```

**参数:**
- `file`: AB1文件（必需）
- `window_size`: 显示窗口大小（可选，默认10bp）
- `ref_file`: 参考基因组文件路径（可选）
- `chrom`: 手动指定染色体（可选）
- `pos`: 手动指定位置（可选）
- `return_json`: 是否返回JSON数据（可选，默认False）

**返回:**
- `return_json=False`: PNG图片文件
- `return_json=True`: JSON数据

**示例（curl）:**
```bash
# 返回图片
curl -X POST "http://localhost:8000/api/align" \
  -F "file=@sample.ab1" \
  -F "window_size=10" \
  -o output.png

# 返回JSON
curl -X POST "http://localhost:8000/api/align" \
  -F "file=@sample.ab1" \
  -F "window_size=10" \
  -F "return_json=true"
```

**示例（Python）:**
```python
import requests

# 返回图片
with open('sample.ab1', 'rb') as f:
    response = requests.post(
        'http://localhost:8000/api/align',
        files={'file': f},
        data={'window_size': 10}
    )
with open('output.png', 'wb') as f:
    f.write(response.content)

# 返回JSON
with open('sample.ab1', 'rb') as f:
    response = requests.post(
        'http://localhost:8000/api/align',
        files={'file': f},
        data={'window_size': 10, 'return_json': True}
    )
result = response.json()
print(result)
```

### 4. 批量文件比对
```
POST /api/align/batch
```

**参数:**
- `files`: AB1文件列表（必需）
- `window_size`: 显示窗口大小（可选，默认10bp）
- `ref_file`: 参考基因组文件路径（可选）
- `return_json`: 是否返回JSON数据（可选，默认False）

**示例（curl）:**
```bash
curl -X POST "http://localhost:8000/api/align/batch" \
  -F "files=@sample1.ab1" \
  -F "files=@sample2.ab1" \
  -F "window_size=10" \
  -F "return_json=true"
```

**示例（Python）:**
```python
import requests

files = [
    ('files', open('sample1.ab1', 'rb')),
    ('files', open('sample2.ab1', 'rb'))
]

response = requests.post(
    'http://localhost:8000/api/align/batch',
    files=files,
    data={'window_size': 10, 'return_json': True}
)

result = response.json()
print(result)
```

### 5. JSON数据返回
```
POST /api/align/json
```

专门用于返回JSON数据的端点，不生成图片。

**参数:**
- `file`: AB1文件（必需）
- `window_size`: 显示窗口大小（可选，默认10bp）
- `ref_file`: 参考基因组文件路径（可选）
- `chrom`: 手动指定染色体（可选）
- `pos`: 手动指定位置（可选）

**返回JSON数据结构:**
```json
{
  "status": "success",
  "filename": "sample.ab1",
  "chrom": "chr7",
  "pos": 140453136,
  "window_size": 10,
  "sample_sequence": "ATCG...",
  "reference_sequence": "ATCG...",
  "variants": [
    {
      "position": 15,
      "type": "mutation",
      "ref": "A",
      "alt": "T"
    }
  ],
  "alignment_score": 98.5,
  "cigar": "100M",
  "mapq": 60,
  "timestamp": "2024-01-01T12:00:00"
}
```

## 返回数据说明

### 变异类型
- `mutation`: 碱基替换
- `deletion`: 缺失
- `insertion`: 插入

### 变异数据结构
```json
{
  "position": 15,
  "type": "mutation",
  "ref": "A",
  "alt": "T"
}
```

## 注意事项

1. 确保已安装BWA工具并配置好参考基因组索引
2. 默认参考基因组路径: `/home/df/test/hg19/hs37d5.fa`
3. 需要同目录下的 `get_reference.py` 脚本
4. 临时文件存储在系统临时目录中

## 依赖项

- FastAPI
- Uvicorn
- Biopython
- NumPy
- Matplotlib
- BWA（系统工具）

## 许可证

MIT License
>>>>>>> 8172053 (Initial commit: Add Sanger sequencing analysis tools)
