#!/usr/bin/env python3
"""
FastAPI应用 - 基因组比对API
提供AB1文件与参考基因组比对的REST API接口
"""

from fastapi import FastAPI, File, UploadFile, Form, HTTPException, Request
from fastapi.responses import JSONResponse, FileResponse, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from typing import Optional, List
import os
import tempfile
import shutil
from datetime import datetime

from plot_genome_alignment import (
    plot_genome_alignment
)

app = FastAPI(
    title="基因组比对API",
    description="AB1文件与参考基因组比对服务",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 临时文件目录
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMP_DIR = os.getenv("TEMP_DIR", tempfile.mkdtemp(prefix="genome_api_"))
OUTPUT_DIR = os.getenv("OUTPUT_DIR", os.path.join(BASE_DIR, "outputs"))
ASSETS_DIR = os.path.join(BASE_DIR, "assets")

os.makedirs(TEMP_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(ASSETS_DIR, exist_ok=True)

# 参考基因组文件路径
GENOME_VERSIONS = {
    "hg19": os.getenv("HG19_REF_FILE", "/data/reference/hg19/hs37d5.fa"),
    "hg38": os.getenv("HG38_REF_FILE", "/data/reference/hg38/hg38.fa")
}

# 验证基因组文件是否存在
AVAILABLE_GENOMES = {}
for version, path in GENOME_VERSIONS.items():
    if os.path.exists(path):
        AVAILABLE_GENOMES[version] = path

app.mount("/outputs", StaticFiles(directory=OUTPUT_DIR), name="outputs")
app.mount("/assets", StaticFiles(directory=ASSETS_DIR), name="assets")


@app.get("/")
async def root():
    """API根路径"""
    return {
        "message": "基因组比对API",
        "version": "1.0.0",
        "available_genomes": list(AVAILABLE_GENOMES.keys()),
        "endpoints": {
            "GET /full_sequence.html": "打开测试页面",
            "POST /api/align": "单个AB1文件比对",
            "POST /api/align/batch": "批量AB1文件比对",
            "GET /api/health": "健康检查"
        }
    }


@app.get("/full_sequence.html")
async def get_full_sequence_page():
    """返回full_sequence.html页面"""
    html_file = os.path.join(os.path.dirname(__file__), "full_sequence.html")
    if not os.path.exists(html_file):
        raise HTTPException(status_code=404, detail="HTML文件未找到")
    return FileResponse(html_file, media_type="text/html")


@app.get("/api/health")
async def health_check():
    """健康检查"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "available_genomes": list(AVAILABLE_GENOMES.keys())
    }


@app.post("/api/align")
async def align_single_file(
    request: Request,
    file: UploadFile = File(...),
    genome_version: str = Form("hg19"),
    window_size: int = Form(10),
    ref_file: Optional[str] = Form(None),
    chrom: Optional[str] = Form(None),
    pos: Optional[int] = Form(None),
    return_json: bool = Form(False)
):
    """
    单个AB1文件比对
    
    参数:
        file: AB1文件
        genome_version: 参考基因组版本（hg19或hg38，默认hg19）
        window_size: 显示窗口大小（默认10bp）
        ref_file: 参考基因组文件路径（可选，优先级高于genome_version）
        chrom: 手动指定染色体（可选）
        pos: 手动指定位置（可选）
        return_json: 是否返回JSON数据（默认False，返回图片）
    """
    try:
        # 确定参考基因组文件
        if ref_file is None:
            if genome_version not in AVAILABLE_GENOMES:
                raise HTTPException(status_code=400, detail=f"不支持的基因组版本: {genome_version}。可用版本: {list(AVAILABLE_GENOMES.keys())}")
            ref_file = AVAILABLE_GENOMES[genome_version]
        
        # 保存上传的文件
        temp_file_path = os.path.join(TEMP_DIR, file.filename)
        with open(temp_file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # 生成输出文件名
        output_filename = f"alignment_{os.path.splitext(file.filename)[0]}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        output_file_path = os.path.join(OUTPUT_DIR, output_filename)
        
        # 执行比对
        result = plot_genome_alignment(
            temp_file_path,
            output_file=output_file_path,
            window_size=window_size,
            ref_file=ref_file,
            chrom=chrom,
            pos=pos
        )
        
        # 清理临时文件
        os.remove(temp_file_path)
        
        if result is None:
            raise HTTPException(status_code=400, detail="比对失败")
        
        if return_json:
            # 返回JSON数据
            return JSONResponse(content={
                "status": "success",
                "filename": file.filename,
                "output_file": output_filename,
                "image_url": str(request.url_for("outputs", path=output_filename)),
                "chrom": result.get("chrom"),
                "pos": result.get("pos"),
                "window_size": window_size,
                "variants": result.get("variants", []),
                "alignment_score": result.get("alignment_score"),
                "timestamp": datetime.now().isoformat()
            })
        else:
            # 返回图片文件
            return FileResponse(
                output_file_path,
                media_type="image/png",
                filename=output_filename
            )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"处理失败: {str(e)}")


@app.post("/api/align/batch")
async def align_batch_files(
    request: Request,
    files: List[UploadFile] = File(...),
    genome_version: str = Form("hg19"),
    window_size: int = Form(10),
    ref_file: Optional[str] = Form(None),
    return_json: bool = Form(True),
    chrom: Optional[str] = Form(None),
    pos: Optional[int] = Form(None)
):
    """
    批量AB1文件比对
    
    参数:
        files: AB1文件列表
        genome_version: 参考基因组版本（hg19或hg38，默认hg19）
        window_size: 显示窗口大小（默认10bp）
        ref_file: 参考基因组文件路径（可选，优先级高于genome_version）
        return_json: 是否返回JSON数据（默认True，返回JSON；False返回ZIP压缩包）
        chrom: 手动指定染色体（可选）
        pos: 手动指定位置（可选）
    """
    try:
        results = []
        output_files = []
        
        for file in files:
            file_result = {
                "filename": file.filename,
                "status": "pending",
                "chrom": None,
                "pos": None,
                "window_size": window_size,
                "sample_sequence": None,
                "reference_sequence": None,
                "variants": [],
                "alignment_score": None,
                "cigar": None,
                "mapq": None,
                "sample_length": None,
                "reference_length": None,
                "error": None
            }
            
            if not file.filename.endswith('.ab1'):
                file_result["status"] = "skipped"
                file_result["error"] = "不是AB1文件"
                results.append(file_result)
                continue
            
            try:
                # 确定参考基因组文件
                current_ref_file = ref_file
                if current_ref_file is None:
                    if genome_version not in AVAILABLE_GENOMES:
                        file_result["status"] = "failed"
                        file_result["error"] = f"不支持的基因组版本: {genome_version}"
                        results.append(file_result)
                        continue
                    current_ref_file = AVAILABLE_GENOMES[genome_version]
                
                # 保存上传的文件
                temp_file_path = os.path.join(TEMP_DIR, file.filename)
                with open(temp_file_path, "wb") as buffer:
                    shutil.copyfileobj(file.file, buffer)
                
                # 生成输出文件名
                output_filename = f"alignment_{os.path.splitext(file.filename)[0]}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
                output_file_path = os.path.join(OUTPUT_DIR, output_filename)
                
                # 执行比对
                result = plot_genome_alignment(
                    temp_file_path,
                    output_file=output_file_path,
                    window_size=window_size,
                    ref_file=current_ref_file,
                    chrom=chrom,
                    pos=pos
                )
                
                # 清理临时文件
                os.remove(temp_file_path)
                
                if result:
                    file_result["status"] = "success"
                    file_result["output_file"] = output_filename
                    file_result["image_url"] = str(request.url_for("outputs", path=output_filename))
                    file_result["chrom"] = result.get("chrom")
                    file_result["pos"] = result.get("pos")
                    file_result["sample_sequence"] = result.get("sample_sequence")
                    file_result["reference_sequence"] = result.get("reference_sequence")
                    file_result["variants"] = result.get("variants", [])
                    file_result["alignment_score"] = result.get("alignment_score")
                    file_result["cigar"] = result.get("cigar")
                    file_result["mapq"] = result.get("mapq")
                    file_result["sample_length"] = result.get("sample_length")
                    file_result["reference_length"] = result.get("reference_length")
                    output_files.append(output_file_path)
                else:
                    file_result["status"] = "failed"
                    file_result["error"] = "比对失败"
                
                results.append(file_result)
                
            except Exception as e:
                file_result["status"] = "failed"
                file_result["error"] = str(e)
                results.append(file_result)
        
        if return_json:
            # 返回JSON数据
            successful_count = len([r for r in results if r["status"] == "success"])
            failed_count = len([r for r in results if r["status"] == "failed"])
            skipped_count = len([r for r in results if r["status"] == "skipped"])
            
            return JSONResponse(content={
                "status": "completed",
                "total_files": len(files),
                "successful": successful_count,
                "failed": failed_count,
                "skipped": skipped_count,
                "results": results,
                "timestamp": datetime.now().isoformat()
            })
        else:
            # 返回ZIP压缩包（包含所有成功的图片）
            import zipfile
            import io
            
            successful_results = [r for r in results if r["status"] == "success"]
            if not successful_results:
                raise HTTPException(status_code=400, detail="所有文件比对失败")
            
            # 创建ZIP文件
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for result in successful_results:
                    output_path = os.path.join(OUTPUT_DIR, result["output_file"])
                    if os.path.exists(output_path):
                        zipf.write(output_path, result["output_file"])
            
            zip_buffer.seek(0)
            
            return Response(
                content=zip_buffer.getvalue(),
                media_type="application/zip",
                headers={
                    "Content-Disposition": f"attachment; filename=alignment_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
                }
            )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"批量处理失败: {str(e)}")


@app.post("/api/align/json")
async def align_to_json(
    file: UploadFile = File(...),
    genome_version: str = Form("hg19"),
    window_size: int = Form(10),
    ref_file: Optional[str] = Form(None),
    chrom: Optional[str] = Form(None),
    pos: Optional[int] = Form(None)
):
    """
    AB1文件比对并返回JSON数据（不生成图片）
    
    参数:
        file: AB1文件
        genome_version: 参考基因组版本（hg19或hg38，默认hg19）
        window_size: 显示窗口大小（默认10bp）
        ref_file: 参考基因组文件路径（可选，优先级高于genome_version）
        chrom: 手动指定染色体（可选）
        pos: 手动指定位置（可选）
    """
    try:
        # 确定参考基因组文件
        if ref_file is None:
            if genome_version not in AVAILABLE_GENOMES:
                raise HTTPException(status_code=400, detail=f"不支持的基因组版本: {genome_version}。可用版本: {list(AVAILABLE_GENOMES.keys())}")
            ref_file = AVAILABLE_GENOMES[genome_version]
        
        # 保存上传的文件
        temp_file_path = os.path.join(TEMP_DIR, file.filename)
        with open(temp_file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # 执行比对（不生成图片）
        result = plot_genome_alignment(
            temp_file_path,
            output_file=None,
            window_size=window_size,
            ref_file=ref_file,
            chrom=chrom,
            pos=pos
        )
        
        # 清理临时文件
        os.remove(temp_file_path)
        
        if result is None:
            raise HTTPException(status_code=400, detail="比对失败")
        
        # 返回详细的JSON数据
        return JSONResponse(content={
            "status": "success",
            "filename": file.filename,
            "chrom": result.get("chrom"),
            "pos": result.get("pos"),
            "window_size": window_size,
            "sample_sequence": result.get("sample_sequence"),
            "reference_sequence": result.get("reference_sequence"),
            "variants": result.get("variants", []),
            "alignment_score": result.get("alignment_score"),
            "cigar": result.get("cigar"),
            "mapq": result.get("mapq"),
            "timestamp": datetime.now().isoformat()
        })

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"处理失败: {str(e)}")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
