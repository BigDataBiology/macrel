package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.Entity.FeedbackInfo;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.utils.EffectiveCheckUtils;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RestController;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 19:47
 * QQ:776748935
 */

@RestController
@CrossOrigin
public class FeedbackController {

    @PostMapping("/feedback")
    public ResultObject feedback(FeedbackInfo feedbackInfo){
        if (!EffectiveCheckUtils.strEffectiveCheck(feedbackInfo.getContent())){
            return ResultObject.failure(ResultCode.FEEDBACK_CONTENT_EMPTY);
        }

        if (EffectiveCheckUtils.strEffectiveCheck(feedbackInfo.getUserEmail())){
            if (!EffectiveCheckUtils.emailEffectiveCheck(feedbackInfo.getUserEmail())){
                return ResultObject.failure(ResultCode.FEEDBACK_EMAIL_FORMAT_ERROR);
            }
        }

        return ResultObject.success(feedbackInfo);
    }
}
