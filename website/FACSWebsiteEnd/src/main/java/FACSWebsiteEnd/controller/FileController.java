package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.EffectiveCheckUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.io.FileSystemResource;
import org.springframework.core.io.InputStreamResource;
import org.springframework.http.HttpHeaders;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import java.io.*;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:55
 * QQ:776748935
 */

@RestController
@CrossOrigin
@RequestMapping("/file")
public class FileController {

    @Autowired
    private FileService fileService;

    @PostMapping("/upload")
    public ResultObject upload(@RequestParam(value = "file") MultipartFile file) {

        if (EffectiveCheckUtils.fileEffectiveCheck(file)){
            return ResultObject.failure(ResultCode.FILE_IS_NULL);
        }

        FileInfo fileInfo = fileService.upload(file);

        if (fileInfo != null){
            return ResultObject.success();
        } else {
            return ResultObject.failure(ResultCode.FILE_SAVE_FAIL);
        }

    }

    @GetMapping("/download")
    public ResponseEntity download() throws IOException {
        String filePath = "H:\\我的图片\\公开\\新浪微博" + "\\img-00d1577eabfa100ac14f643187343e3b.jpg";
        FileSystemResource file = new FileSystemResource(filePath);
        HttpHeaders headers = new HttpHeaders();
        // 设置下载文件默认的名称
        headers.add("Content-Disposition","attachment;filename=123.jpg");

        return ResponseEntity.ok()
                .headers(headers)
                .contentLength(file.contentLength())
                // 表示响应的内容是通过字节流的方式进行传输的
                .contentType(MediaType.parseMediaType("application/octet-stream"))
                .body(new InputStreamResource(file.getInputStream()));
    }
}
