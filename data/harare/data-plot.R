d = read.csv("harareClean.csv")
d$date = as.Date(d$date,"%m/%d/%Y")
d$new = diff(c(0,d$total))

library(ggplot2)
qplot(date, new, data=d, geom="point")
qplot(date, total, data=d, geom="point")


meeting  = as.Date("01-05-2010","%d-%m-%Y")
campaign = as.Date(c("24-05-2010","02-06-2010"),"%d-%m-%Y")


p = ggplot(d, aes(x = date, y = new))+geom_point()
q = ggplot(d, aes(x = date, y = total)) + geom_point()
meeting.line = geom_vline(xintercept = as.numeric(meeting), col="red")
campaign.rect = geom_rect(xmin=as.numeric(campaign[1]), xmax=as.numeric(campaign[2]), ymin=-100, ymax=1000, fill="green")

pdf("new-cases.pdf")
print(p+meeting.line+campaign.rect+geom_point()+ labs(title = "New weekly cases", y="Number of new cases", x="Date"))
dev.off()

pdf("cumulative-cases.pdf")
print(q+meeting.line+campaign.rect+geom_point()+ labs(title = "Cumulative cases", y="Cumulative number of cases", x="Date"))
dev.off()


q(ifelse(interactive(),"ask","no"))

