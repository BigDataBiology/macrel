package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.DataUploaded;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.utils.EffectiveCheckUtils;
import FACSWebsiteEnd.utils.FileUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.multipart.MultipartFile;

import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/11/29 11:09
 * QQ:776748935
 */

@RestController
@CrossOrigin
@RequestMapping("/facs")
public class FacsController {

    @Autowired
    private FacsService facsService;

    @PostMapping("/prediction")
    public ResultObject analysis(DataUploaded dataUploaded){

        FileInfo fileInfo;

        // 校验上传的文本和文件
        if (!EffectiveCheckUtils.strEffectiveCheck(dataUploaded.getSequence())
                && !EffectiveCheckUtils.fileEffectiveCheck(dataUploaded.getFile())){
            return ResultObject.failure(ResultCode.DATA_IS_EMPTY);
        }

        // 上传的是文本
        if (EffectiveCheckUtils.strEffectiveCheck(dataUploaded.getSequence())){
            fileInfo = facsService.saveSequenceToFile(dataUploaded.getSequence());
        } else {
            // 上传的是文件
            MultipartFile file = dataUploaded.getFile();
            Map fileInformation = FileUtils.getFileInformation(file);

            String type = fileInformation.get("type").toString();

            if (type != null){
                // 判断是否是指定类型的文件
                if (Constant.EXTENSION.equals(type)){
                    fileInfo = facsService.saveFile(file);
                }else {
                    return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTQ_ERROR);
                }
            } else {
                return ResultObject.failure(ResultCode.FILETYPE_UNKNOWN_ERROR);
            }

        }

        //todo 调用脚本
        String sequenceType = "protein";
        String mode = "r";
        String read_1 = "pathOfRead_1";
        facsService.callShell(sequenceType,mode,read_1);

        return ResultObject.success(fileInfo);
    }

}
