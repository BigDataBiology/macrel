package FACSWebsiteEnd.common;

/**
 * @author HiramHe
 */
public class ResultObject {

    private Integer code;
    private String msg;
    private Object data;

    public static final String NO_DATA = "NO_DATA";

    public ResultObject() {}

    public ResultObject(Integer code, String msg) {
        this.code = code;
        this.msg = msg;
    }

    public static ResultObject success() {
        ResultObject result = new ResultObject();
        result.setResultCode(ResultCode.SUCCESS);
        result.setData(NO_DATA);
        return result;
    }

    public static ResultObject success(Object data) {
        ResultObject result = new ResultObject();
        result.setResultCode(ResultCode.SUCCESS);
        result.setData(data);
        return result;
    }

    public static ResultObject failure(ResultCode resultCode) {
        ResultObject result = new ResultObject();
        result.setResultCode(resultCode);
        return result;
    }

    public static ResultObject failure(ResultCode resultCode, Object data) {
        ResultObject result = new ResultObject();
        result.setResultCode(resultCode);
        result.setData(data);
        return result;
    }

    public void setResultCode(ResultCode code) {
        this.code = code.code();
        this.msg = code.message();
    }

    public Integer getCode() {
        return code;
    }

    public void setCode(Integer code) {
        this.code = code;
    }

    public String getMsg() {
        return msg;
    }

    public void setMsg(String msg) {
        this.msg = msg;
    }

    public Object getData() {
        return data;
    }

    public void setData(Object data) {
        this.data = data;
    }
}
