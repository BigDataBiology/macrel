package FACSWebsiteEnd.controller;

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
    private PipelineProperties pipelineConfiguration;
    @Autowired
    private RemoteProperties remoteConfiguration;

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

        // 判断website与pipeline是否在同一台服务器上
        if (!remoteConfiguration.getEnableRemote()){
            FacsUtils.createFolderLocally(pipelineConfiguration.getInputDir());
            FacsUtils.createFolderLocally(pipelineConfiguration.getOutputDir());
        } else {
//            Session session = CommandUtils.connect(remoteConfiguration);
            FacsUtils.createFolderRemotely(remoteConfiguration,pipelineConfiguration);
        }

        FileInfo fileInfo = null;
        // 保存数据
        // 上传的是文本
        if (EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())){
            if (!remoteConfiguration.getEnableRemote()){
                fileInfo = fileService.saveTextToFileLocally(predictionForm.getTextData(),pipelineConfiguration.getInputDir(),Constant.FA);
            } else {
                //todo
            }
        } else {
            // 上传的是文件
            MultipartFile file = predictionForm.getFile();
            Map fileInformation = FileUtils.getFileInformation(file);

            String extension = fileInformation.get("extension").toString();
            if (extension != null){
                // 判断是否是指定类型的文件,格式需为 fasta、fa
                if (Constant.FASTA.equals(extension) || Constant.FA.equals(extension)){
                    if (!remoteConfiguration.getEnableRemote()){
                        fileInfo = fileService.uploadFileToLocal(file,pipelineConfiguration.getInputDir());
                    } else {
                        //todo
                    }
                } else {
                    return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTA_OR_FA_ERROR);
                }
            } else {
                return ResultObject.failure(ResultCode.FILETYPE_UNKNOWN_ERROR);
            }
        }
//        System.out.println(fileInfo);

        // 在Linux上把当前输出的文件夹创建出来
        String folderName = null;
        if (fileInfo != null){
            folderName = fileInfo.getFilenameWithOutExtension();
        }
        String currentOutputDir = null;
        if (!remoteConfiguration.getEnableRemote()){
            currentOutputDir = FacsUtils.createFolderLocally(pipelineConfiguration.getOutputDir(),folderName);
        } else {
            //todo
        }

        // just for test
//        String inputFilePath = "/home/HiramHe/facs_data_uploadByUser/sequence-7ead845137a64d08b0092b8224766e25.fa";
//        fileInfo.setPath(inputFilePath);
//        String currentOutDir = "/home/HiramHe/facs_out/sequence-7ead845137a64d08b0092b8224766e25";

        // 调用pipeline，对数据进行处理
        String outputFilePath = null;
        if (fileInfo!= null){
            outputFilePath = facsService.callPipeline(pipelineConfiguration.getHome(), fileInfo, currentOutputDir, dataType, remoteConfiguration.getEnableRemote());
        }

        // 读取结果
        List<Object> objects = null;
        if (!remoteConfiguration.getEnableRemote()){
            objects = facsService.readResultsLocally(outputFilePath);
        } else {
            //todo
        }

        // 封装数据
        PredictionOut predictionOut = new PredictionOut();
        predictionOut.setObjects(objects);
        predictionOut.setFilePath(outputFilePath);

        // 返回数据
        return ResultObject.success(predictionOut);

    }

}
