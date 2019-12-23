package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.Entity.PredictionForm;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.Entity.PredictionOut;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
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
@CrossOrigin(origins = "*")
@RequestMapping("/facs")
public class FacsController {

    @Autowired
    private FacsService facsService;
    @Autowired
    private FileService fileService;

    @Autowired
    private PipelineProperties pipelineProperties;
    @Autowired
    private RemoteProperties remoteProperties;

    @PostMapping("/prediction")
    public ResultObject analysis(PredictionForm predictionForm){

        // just for test
//        String savedDir = Constant.FILESAVED_WIN_DIR;

        String dataType = predictionForm.getDataType();

        // 校验数据类型是否为空
        if (!EffectiveCheckUtils.strEffectiveCheck(dataType)){
            return ResultObject.failure(ResultCode.DATATYPE_EMPTY);
        }

        // 校验上传的文本和文件
        if (!EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())
                && !EffectiveCheckUtils.fileEffectiveCheck(predictionForm.getFile())){
            return ResultObject.failure(ResultCode.DATA_IS_EMPTY);
        }

        FileInfo fileInfo = null;
        // 保存数据
        // 上传的是文本
        if (EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())){
            String content = predictionForm.getTextData();
            if (!remoteProperties.getEnableRemote()){
                fileInfo = fileService.saveTextToFileLocally(content,pipelineProperties.getInputDir(),Constant.FA);
            } else {
                // done
                fileInfo = fileService.saveContentToFileRemotely(remoteProperties, pipelineProperties.getInputDir(), Constant.FA, content);
            }
        } else {
            // 上传的是文件
            MultipartFile file = predictionForm.getFile();
            Map fileInformation = FileUtils.getFileInformation(file);

            Object extension = fileInformation.get("extension");
            if (extension != null){
                extension = extension.toString();
                // 判断是否是指定类型的文件,格式需为 fasta、fa
                if (Constant.FASTA.equals(extension) || Constant.FA.equals(extension)){
                    if (!remoteProperties.getEnableRemote()){
                        fileInfo = fileService.uploadFileToLocal(file,pipelineProperties.getInputDir());
                    } else {
                        fileInfo = fileService.uploadFileToRemote(remoteProperties, pipelineProperties.getInputDir(), file);
                    }
                } else {
                    return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTA_OR_FA_ERROR);
                }
            } else {
                return ResultObject.failure(ResultCode.FILETYPE_UNKNOWN_ERROR);
            }
        }
//        System.out.println(fileInfo);

        if (fileInfo == null){
            return ResultObject.failure(ResultCode.FILE_INFO_FAIL);
        }

        // 在Linux上把当前的输出文件夹创建出来
        String folderName = fileInfo.getFilenameWithOutExtension();
        String currentOutputDir = null;
        if (!remoteProperties.getEnableRemote()){
            currentOutputDir = FacsUtils.createFolderLocally(pipelineProperties.getOutputDir(),folderName);
        } else {
            // done
            currentOutputDir = FacsUtils.createFolderRemotely(remoteProperties,pipelineProperties.getOutputDir(),folderName);
        }
//        System.out.println("currentOutputDir:"+currentOutputDir);

        if (currentOutputDir == null){
            return ResultObject.failure(ResultCode.FOLDER_CREATE_FAIL);
        }

        // just for test
//        String inputFilePath = "/home/HiramHe/facs_data_uploadByUser/sequence-7ead845137a64d08b0092b8224766e25.fa";
//        fileInfo.setPath(inputFilePath);
//        String currentOutDir = "/home/HiramHe/facs_out/sequence-7ead845137a64d08b0092b8224766e25";

        // 调用pipeline，对数据进行处理
        String outputFilePath = facsService.callPipeline(pipelineProperties.getHome(), fileInfo, currentOutputDir, dataType, remoteProperties);
//        System.out.println("outputFilePath:"+outputFilePath);

        if (outputFilePath == null){
            return ResultObject.failure(ResultCode.PIPELINE_CALL_ERROR);
        }

        // 读取结果
        List<Object> objects = null;
        if (!remoteProperties.getEnableRemote()){
            objects = fileService.readFileToObjectFromLocal(outputFilePath,new FacsOutTsv());
        } else {
            // done
            objects = fileService.readFileToObjectFromRemote(remoteProperties,outputFilePath,new FacsOutTsv());
        }

        if (objects == null){
            return ResultObject.failure(ResultCode.RESULT_READ_FAIL);
        }

        // 封装数据
        PredictionOut predictionOut = new PredictionOut();
        predictionOut.setObjects(objects);
        predictionOut.setFilePath(outputFilePath);

        // 返回数据
        return ResultObject.success(predictionOut);

    }

}
