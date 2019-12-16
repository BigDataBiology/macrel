package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.PredictionForm;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.Entity.PredictionOut;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.system.ApplicationHome;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.multipart.MultipartFile;

import java.io.File;
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

    private ApplicationHome applicationHome = new ApplicationHome(getClass());
    private File jarFile = applicationHome.getSource();
    private String jarPath = jarFile.getParentFile().toString();

    @PostMapping("/prediction")
    public ResultObject analysis(PredictionForm predictionForm){

        String savedDir = FacsUtils.makeSavedFolderOnLinux(jarPath);
        String allOutDir = FacsUtils.makeAllOutFolderOnLinux(jarPath);
        // todo:just for test
//        String savedDir = Constant.FILESAVED_WIN_DIR;

        FileInfo fileInfo = null;
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

        // 保存数据
        // 上传的是文本
        if (EffectiveCheckUtils.strEffectiveCheck(predictionForm.getTextData())){
            fileInfo = fileService.saveTextToFile(predictionForm.getTextData(),savedDir,Constant.FA);
        } else {
            // 上传的是文件
            MultipartFile file = predictionForm.getFile();
            Map fileInformation = FileUtils.getFileInformation(file);

            String extension = fileInformation.get("extension").toString();

            if (extension != null){
                // 判断是否是指定类型的文件,格式需为 fasta、fa

                if (Constant.FASTA.equals(extension) || Constant.FA.equals(extension)){
                    fileInfo = fileService.uploadFileToLocal(file,savedDir);
                } else {
                    return ResultObject.failure(ResultCode.FILETYPE_NOT_FASTA_OR_FA_ERROR);
                }
            } else {
                return ResultObject.failure(ResultCode.FILETYPE_UNKNOWN_ERROR);
            }

        }
//        System.out.println(fileInfo);

        // 调用脚本
        // windows上远程测试时用,在Linux上运行时请注释掉
//        fileInfo.setFullpath(Constant.FILESAVED_REMOTE_DIR + fileInfo.getFilename());

        // 在Linux上把当前输出的文件夹创建出来
        String outfolderName = fileInfo.getFilenameWithOutExtension();
        String currentOutDir = FacsUtils.makeCurrentOutFolderOnLinux(allOutDir,outfolderName);

        // 调用pipeline，对数据进行处理
        // todo: just for test
//        String inputFilePath = "/home/HiramHe/facs_data_uploadByUser/sequence-7ead845137a64d08b0092b8224766e25.fa";
//        fileInfo.setPath(inputFilePath);
//        String currentOutDir = "/home/HiramHe/facs_out/sequence-7ead845137a64d08b0092b8224766e25";

        facsService.callShellScript(fileInfo, currentOutDir, dataType,false);

        String outputFilePath = currentOutDir + Constant.FACS_OUT_FILENAME;
        // 读取结果
        List<Object> objects = facsService.readLocalResults(outputFilePath);

        // 封装数据
        PredictionOut predictionOut = new PredictionOut();
        predictionOut.setObjects(objects);
        predictionOut.setFilePath(outputFilePath);

        // 返回数据
        return ResultObject.success(predictionOut);
    }

}
