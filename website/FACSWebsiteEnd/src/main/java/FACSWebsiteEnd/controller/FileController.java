package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.EffectiveCheckUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;

import javax.servlet.http.HttpServletResponse;
import java.io.*;
import java.net.URLEncoder;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:55
 * QQ:776748935
 */

@RestController
@CrossOrigin(origins = "*")
@RequestMapping("/file")
public class FileController {

    @Autowired
    private FileService fileService;

    @PostMapping("/upload")
    public ResultObject upload(@RequestParam(value = "file") MultipartFile file) {

        if (!EffectiveCheckUtils.fileEffectiveCheck(file)){
            System.out.println("no file.");
            return ResultObject.failure(ResultCode.FILE_IS_NULL);
        }

        FileInfo fileInfo = fileService.uploadFileToLocal(file, Constant.FILESAVED_WIN_DIR);

        if (fileInfo != null){
            return ResultObject.success();
        } else {
            return ResultObject.failure(ResultCode.FILE_SAVE_FAIL);
        }

    }

    @GetMapping("/download")
    public ResultObject download(@RequestParam("filePath") String filePath, HttpServletResponse response) {

        if (!EffectiveCheckUtils.strEffectiveCheck(filePath)){
            return ResultObject.failure(ResultCode.FILE_NOT_EXIST);
        }

        String filenameWithExtension = filePath.lastIndexOf("/") != -1 ?
                filePath.substring(filePath.lastIndexOf("/")+1) : filePath;

        File file = new File(filePath);

        if (file.exists()){
            // 实现文件下载

            byte[] buffer = new byte[1024];
            FileInputStream fileInputStream = null;
            BufferedInputStream bufferedInputStream = null;
            OutputStream outputStream = null;

            try {

                //1.设置文件ContentType类型，这样设置，会自动判断下载文件类型
                response.setContentType("multipart/form-data");
                //2.设置文件头：最后一个参数是设置下载文件名
                String fileName = URLEncoder.encode(filenameWithExtension, "UTF-8");
                response.setHeader("Content-Disposition", "attachment;filename=" + fileName);
                // 便于前端获取文件名
                response.setHeader("fileName", fileName);
                response.setHeader("Access-Control-Expose-Headers", "fileName");
                fileInputStream = new FileInputStream(file);
                bufferedInputStream = new BufferedInputStream(fileInputStream);

                //3.通过response获取OutputStream对象
                outputStream = response.getOutputStream();

                int i = bufferedInputStream.read(buffer);
                while (i != -1) {
                    outputStream.write(buffer, 0, i);
                    i = bufferedInputStream.read(buffer);
                }

                return null;

            } catch (IOException e) {
                return null;
            } finally {
                if (bufferedInputStream != null){
                    try {
                        bufferedInputStream.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                if (fileInputStream != null){
                    try {
                        fileInputStream.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }

        } else {
            return ResultObject.failure(ResultCode.FILE_NOT_EXIST);
        }

    }
}
