#include "clang/AST/ASTContext.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendAction.h"
#include "clang/Rewrite/Core/Rewriter.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "llvm/Support/CommandLine.h"

using namespace clang;
using namespace clang::tooling;
static llvm::cl::OptionCategory MyToolCategory("my-tool options");

using namespace clang::ast_matchers;

class TestHandler : public MatchFinder::MatchCallback {
public:
    TestHandler(Rewriter &R) : TheRewriter(R) {};
    virtual void run(const MatchFinder::MatchResult &Result) {
        ASTContext *Context = Result.Context;
        llvm::outs() << "Got here\n";
        const auto ActorExpr = Result.Nodes.getNodeAs<CXXMemberCallExpr>("call");
        ActorExpr->dump();
    }
private:
    Rewriter &TheRewriter;
};

class SelectorLambdaHandler : public MatchFinder::MatchCallback {
public:
    SelectorLambdaHandler(Rewriter &R) : TheRewriter(R) {};

    virtual void run(const MatchFinder::MatchResult &Result) {
        ASTContext *Context = Result.Context;
        // lambda
        auto const *Lambda = Result.Nodes.getNodeAs<CXXRecordDecl>("Lambda");
        auto const *mDecl = Lambda->getLambdaCallOperator();
        // hs_ptr->send();
        const auto ActorSendExpr = Result.Nodes.getNodeAs<CXXMemberCallExpr>("send");

        // Analyze the captured vars by the lambda
        // 1. Non pointer type: treated as an index to global arrays
        // 2. Pointer type: treated as a global array
        // TODO: maybe some selection is required
        std::string packet_type;
        std::vector<const VarDecl *> nonptrvars;
        std::vector<const VarDecl *> ptrvars;
        if (mDecl->hasBody()) {
            llvm::outs() << "[Current Lambda]\n";
            mDecl->getBody()->dump();
            llvm::outs() << Lexer::getSourceText(CharSourceRange::getTokenRange(mDecl->getBody()->getSourceRange()), TheRewriter.getSourceMgr(), TheRewriter.getLangOpts());
            for (auto const &Capture : Lambda->captures()) {
                auto const var = Capture.getCapturedVar();
                if (!var->getType().getTypePtr()->isPointerType()) {
                    nonptrvars.push_back(var);
                } else {
                    ptrvars.push_back(var);
                }
            }
        }

        if (nonptrvars.size() > 1) {
            llvm::outs() << "Warning: # of non-ptr vars is > 1\n";
        }
        packet_type = nonptrvars[0]->getType().getAsString();

        // Synthesize a selector class
        const DeclRefExpr *ActorInstantiation = nullptr;
        for (auto const &p : ActorSendExpr->getImplicitObjectArgument()->children()) {
            ActorInstantiation = dyn_cast<DeclRefExpr>(p);
            if (ActorInstantiation) {
                break;
            }
        }
        if (ActorInstantiation) {
            ActorInstantiation->getDecl()->dump();
            const VarDecl *ActorInstance = dyn_cast<VarDecl>(ActorInstantiation->getDecl());
            // Synthesize
            {
                std::string mes = "\n";
                mes += "class SynthesizedActor : public hclib::Actor<" + packet_type + "> { \n";
                mes += "public: \n";
                for (auto const v : ptrvars) {
                    mes += v->getType().getAsString() + " " + v->getName().str() + ";\n";
                }
                mes += "void process(";
                for (auto const v: nonptrvars) {
                    mes += v->getType().getAsString() + " " + v->getName().str();
                    mes += ", ";
                }
                mes += "int sender_rank)";
                mes += Lexer::getSourceText(CharSourceRange::getTokenRange(mDecl->getBody()->getSourceRange()), TheRewriter.getSourceMgr(), TheRewriter.getLangOpts());
                mes += "\n";
                mes += "SynthesizedActor(";
                for (auto const v : ptrvars) {
                    mes += v->getType().getAsString() + " _" + v->getName().str();
                    if (v != *(ptrvars.end() - 1)) mes += ", ";
                }
                mes += "): ";
                for (auto const v : ptrvars) {
                    mes += v->getName().str() + "(_" + v->getName().str() + ")";
                    if (v != *(ptrvars.end() - 1)) mes += ", ";
                    else mes += " ";
                }
                // mb[0].process = [this](int64_t pkt, int sender_rank) { this->process(pkt, sender_rank);};
                mes += "{\n";
                // TODO: support multiple mailboxes
                mes += "mb[0].process = [this](" + packet_type + " pkt, int sender_rank) { this->process(pkt, sender_rank); };\n";
                mes += "}\n";
                mes += "};\n";
                mes += "SynthesizedActor *" + ActorInstance->getName().str();
                mes += " = new SynthesizedActor(";
                for (auto const v : ptrvars) {
                    mes += v->getName().str();
                    if (v != *(ptrvars.end() - 1)) mes += ", ";
                }
                mes += ");\n";
                mes += "//";
                //ActorInstantiation->getDecl()->getBeginLoc().print(llvm::outs(), TheRewriter.getSourceMgr());
                TheRewriter.InsertText(ActorInstantiation->getDecl()->getBeginLoc(), mes, true, true);
            }
            // Update actor->send
            {
                // TODO: the arguments to send is always int + lambda?
                const DeclRefExpr *arg0 = nullptr;
                for (auto const &p : ActorSendExpr->getArg(0)->children()) {
                    for (auto const &q : p->children()) {
                        arg0 = dyn_cast<DeclRefExpr>(q);
                        if (arg0) {
                            llvm::outs() << "hoge2\n";
                            arg0->dump();
                        }
                    }
                }
                if (arg0) {
                    //const VarDecl argd = dyn_cast<VarDecl>(arg);
                    std::string mes = ActorInstance->getName().str() + "->send(";
                    // col
                    for (auto const v: nonptrvars) {
                        mes += v->getName().str();
                        mes += ", ";
                    }
                    // pe
                    mes += arg0->getDecl()->getName().str() + ");\n ";
                    mes += "//";
                    TheRewriter.InsertText(ActorSendExpr->getBeginLoc(), mes, true, true);
                }
            }
        } else {
            llvm::outs() << "ActorInstatiation not found\n";
        }
    }
private:
    Rewriter &TheRewriter;
};

class SelectorTransConsumer : public clang::ASTConsumer {
public:
    explicit SelectorTransConsumer(Rewriter &R)
        : HandlerForLambda(R), HandlerForTest(R) {

        auto Lambda = expr(hasType(cxxRecordDecl(isLambda()).bind("Lambda")));
        Matcher.addMatcher(cxxMemberCallExpr(callee(cxxMethodDecl(hasName("send"))),
                                             hasArgument(1, Lambda)).bind("send"),
                           &HandlerForLambda);
    }

    virtual void HandleTranslationUnit(clang::ASTContext &Context) {
        // Run the matchers when we have the whole TU parsed.
        Matcher.matchAST(Context);
    }
private:
    SelectorLambdaHandler HandlerForLambda;
    TestHandler HandlerForTest;
    MatchFinder Matcher;
};

class SelectorTransAction : public clang::ASTFrontendAction {
public:
    void EndSourceFileAction() override {
        SourceManager &SM = TheRewriter.getSourceMgr();
        llvm::errs() << "** EndSourceFileAction for: "
                     << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";
        TheRewriter.getEditBuffer(SM.getMainFileID()).write(llvm::outs());
    }

    virtual std::unique_ptr<clang::ASTConsumer> CreateASTConsumer(
        clang::CompilerInstance &Compiler, llvm::StringRef InFile) {
        llvm::errs() << "** Creating AST consumer for: " << InFile << "\n";
        TheRewriter.setSourceMgr(Compiler.getSourceManager(), Compiler.getLangOpts());
        return std::unique_ptr<clang::ASTConsumer>(
            new SelectorTransConsumer(TheRewriter));
    }
private:
    Rewriter TheRewriter;
};

int main(int argc, const char **argv) {
    CommonOptionsParser OptionsParser(argc, argv, MyToolCategory);
    ClangTool Tool(OptionsParser.getCompilations(),
                   OptionsParser.getSourcePathList());

    return Tool.run(newFrontendActionFactory<SelectorTransAction>().get());
}
