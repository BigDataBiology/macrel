package FACSWebsiteEnd.common;

import java.util.ArrayList;
import java.util.List;

/**
 * @author HiramHe
 */
public enum ResultCode {

    /* 成功状态码 */
    SUCCESS(1,"成功"),

    /* 服务器相关 */
    SERVER_INTERNAL_ERROR(500,"服务器内部错误."),

    /* prediction相关：10001-19999 */
    FILE_SAVE_SUCCESS(10001, "文件上传成功."),
    FILE_IS_NULL(10002, "没有选择文件或文件内容为空."),
    FILE_SAVE_FAIL(10003, "内部错误,文件保存失败."),
    DATA_IS_EMPTY(10004, "至少输入序列文本或者上传文件"),
    FILETYPE_NOT_FASTQ_ERROR(10005, "非fastq文件."),
    FILETYPE_UNKNOWN_ERROR(10006, "未知文件."),
    DATATYPE_UNSPECIFIED(10007, "未指定数据类型."),
    FILETYPE_NOT_FASTA_OR_FA_ERROR(10008, "上传的文件格式需为fasta或者fa."),
    FILETYPE_ERROR(10009, "文件类型错误."),
    DATATYPE_EMPTY(10010, "请选择数据类型."),
    FILE_NOT_EXIST(10011, "文件不存在."),

    /* 反馈相关： 20001-29999*/
    FEEDBACK_SUCCESS(20001,"反馈成功"),
    FEEDBACK_CONTENT_EMPTY(20002,"反馈内容为空"),
    FEEDBACK_EMAIL_FORMAT_ERROR(2003,"邮箱格式错误");

    private Integer code;
    private String message;

    ResultCode() {
    }

    ResultCode(Integer code, String message) {
        this.code = code;
        this.message = message;
    }

    public Integer code() {
        return this.code;
    }

    public String message() {
        return this.message;
    }

    public static String getMessage(String name) {
        for (ResultCode item : ResultCode.values()) {
            if (item.name().equals(name)) {
                return item.message;
            }
        }
        return name;
    }

    public static Integer getCode(String name) {
        for (ResultCode item : ResultCode.values()) {
            if (item.name().equals(name)) {
                return item.code;
            }
        }
        return null;
    }

    @Override
    public String toString() {
        return this.name();
    }

    //校验重复的code值
    public static void main(String[] args) {
        ResultCode[] apiResultCodes = ResultCode.values();
        List<Integer> codeList = new ArrayList<>();
        for (ResultCode apiResultCode : apiResultCodes) {
            if (codeList.contains(apiResultCode.code)) {
                System.out.println(apiResultCode.code);
            } else {
                codeList.add(apiResultCode.code());
            }
        }
    }
}
