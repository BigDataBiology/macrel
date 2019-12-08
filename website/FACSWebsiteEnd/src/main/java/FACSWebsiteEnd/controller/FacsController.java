package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.PredictionForm;
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

import java.util.List;
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
    public ResultObject analysis(PredictionForm predictionForm){

        FileInfo fileInfo = null;
        String dataType = predictionForm.getDataType();

        if (!EffectiveCheckUtils.strEffectiveCheck(dataType)){
            return ResultObject.failure(ResultCode.DATATYPE__EMPTY);
        }

        // 校验上传的文本和文件
        if (!EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())
                && !EffectiveCheckUtils.fileEffectiveCheck(predictionForm.getFile())){
            return ResultObject.failure(ResultCode.DATA_IS_EMPTY);
        }

        // 上传的是文本
        if (EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())){
            fileInfo = facsService.saveSequenceToFile(predictionForm.getTextData(),predictionForm.getDataType());
        } else {
            // 上传的是文件
            MultipartFile file = predictionForm.getFile();
            Map fileInformation = FileUtils.getFileInformation(file);

            String extension = fileInformation.get("extension").toString();

            if (extension != null){
                // 判断是否是指定类型的文件

                // PEPTIDES，上传文件的格式需为 fasta、fa
                if (Constant.PEPTIDES.equals(dataType)){
                    if (Constant.FASTA.equals(extension) || Constant.FA.equals(extension)){
                        fileInfo = facsService.saveFile(file);
                    }else {
                        return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTA_OR_FA_ERROR);
                    }
                } else if (Constant.NUCLEOTIDE.equals(dataType)){
                    // NUCLEOTIDE，上传的文件格式需为 fastaq
                    if (Constant.FASTTQ.equals(extension)){
                        fileInfo = facsService.saveFile(file);
                    }else {
                        return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTQ_ERROR);
                    }
                } else {
                    return ResultObject.failure(ResultCode.FILETYPE__ERROR);
                }


            } else {
                return ResultObject.failure(ResultCode.FILETYPE_UNKNOWN_ERROR);
            }

        }

        // 校验数据类型是否为空
        if (!EffectiveCheckUtils.strEffectiveCheck(predictionForm.getDataType())){
            return ResultObject.failure(ResultCode.DATATYPE_UNSPECIFIED);
        }

        List<Object> objects = facsService.callShell(fileInfo, dataType);

        return ResultObject.success(objects);
    }

}
