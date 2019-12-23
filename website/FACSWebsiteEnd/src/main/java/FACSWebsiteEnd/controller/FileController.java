package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;
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
    @Autowired
    RemoteProperties remoteProperties;
    @Autowired
    PipelineProperties pipelineProperties;


    @PostMapping("/upload")
    public ResultObject upload(@RequestParam(value = "file") MultipartFile file) {

        if (!EffectiveCheckUtils.fileEffectiveCheck(file)){
            System.out.println("no file.");
            return ResultObject.failure(ResultCode.FILE_IS_NULL);
        }

        String dir = Constant.FILESAVED_WIN_DIR;
        FileInfo fileInfo = fileService.uploadFileToLocal(file,dir);

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

        ResultObject resultObject = null;
        if (!remoteProperties.getEnableRemote()){
            if (file.exists()){
                resultObject = fileService.getFileForDownloadFromLocal(file,filenameWithExtension,response);
                return resultObject;
            } else {
                return ResultObject.failure(ResultCode.FILE_NOT_EXIST);
            }

        } else {
            // done 远程下载
            resultObject = fileService.getFileForDownloadFromRemote(remoteProperties,filePath,filenameWithExtension,response);
            return resultObject;
        }

    }
}
