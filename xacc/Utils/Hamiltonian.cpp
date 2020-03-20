#include "Hamiltonian.hpp"
#include <algorithm>
#include <cctype>
#include <iostream>
#include "xacc_service.hpp"
#include "expression_parsing_util.hpp"
#include <cassert>
#include "xacc.hpp"
#include "json.hpp"
#include <regex>
#include "Functor.hpp"
#include "Executor.hpp"

namespace {
    std::string removeWhiteSpaces(const std::string& in_str)
    {
        auto result = in_str;
        result.erase(std::remove_if(result.begin(), result.end(), [](unsigned char c){ return std::isspace(c); }), result.end());

        return result;
    }

    std::string toUpperCase(const std::string& in_str)
    {
        auto result = in_str;
        transform(result.begin(), result.end(), result.begin(), ::toupper); 
        return result;
    }

    bool isNumberString(const std::string& in_str)
    {
        return std::find_if(in_str.begin(), in_str.end(), [](const auto in_char){ return in_char != '0' && in_char != '1' &&
                                                                                         in_char != '2' && in_char != '3' &&
                                                                                         in_char != '4' && in_char != '5' &&
                                                                                         in_char != '6' && in_char != '7' &&
                                                                                         in_char != '8' && in_char != '9'; }) == in_str.end();
    }

    bool GetLastOperator(const std::string& in_str, QuaC::QubitOp& out_Op, std::string& out_subString)
    {
        // Find the last '*' character
        const auto pos = in_str.find_last_of("*");
        if (pos == std::string::npos)
        {
            return false;
        }

        const auto opName = toUpperCase(in_str.substr(pos + 1));     
        const auto qubitIdxIter = std::find_if(opName.begin(), opName.end(), [](const auto in_char){ return !isalpha(in_char); } );
        
        if (qubitIdxIter == opName.end())
        {
            return false;
        }

        const auto opStr = opName.substr(0, std::distance(opName.begin(), qubitIdxIter));


        if (opStr.length() < 1 ||  QuaC::ConvertOperatorFromString(opStr) == QuaC::Operator::NA)
        {
            return false;
        }

        const auto qubitIdxStr = opName.substr(std::distance(opName.begin(), qubitIdxIter));
        if (!isNumberString(qubitIdxStr))
        {
            return false;
        }

        const auto op = QuaC::ConvertOperatorFromString(opStr);
        const size_t qIdx = std::stoi(qubitIdxStr);
        
        out_Op = std::make_pair(op, qIdx);
        out_subString = in_str.substr(0, pos);
        return true;
    }
    // Unwrap expression of type blabla*(A +/- B) into blabla*A +/- blabla*B 
    std::vector<std::string> UnwrapOpExpresion(const std::string& in_str)
    {
        if (in_str.back() != ')')
        {
            return {};
        }

        const auto pos = in_str.find_last_of("(");
        
        if (pos == std::string::npos)
        {
            return {};
        }

        std::string coeffExpr = in_str.substr(0, pos);
        std::string wrappedExpr = in_str.substr(pos + 1);
        // remove the last ')' character
        wrappedExpr.pop_back();

        const auto check1 = wrappedExpr.find_first_of("(");
        const auto check2 = wrappedExpr.find_first_of(")");
        if (check1 != std::string::npos || check2 != std::string::npos)
        {
            // Perhaps a nested one, we cannot parse atm
            return {};
        }
        const auto pmPos = wrappedExpr.find_first_of("+") == std::string::npos ? wrappedExpr.find_first_of("-") : wrappedExpr.find_first_of("+");
        
        if (pmPos != std::string::npos)
        {
            std::string expr1 = wrappedExpr.substr(0, pmPos);
            std::string expr2 = wrappedExpr.substr(pmPos + 1);
            // (1.0)* or (-1.0)*
            const std::string signExpr = "(" + std::string(1, wrappedExpr.at(pmPos)) + "1.0)*";
            return { coeffExpr + expr1, signExpr + coeffExpr  + expr2 };
        }

        return {};
    }

    bool TryEvaluateExpression(const std::string& in_exprString, const QuaC::VarsMap& in_vars, double& out_result)
    {
        auto parsingUtil = xacc::getService<xacc::ExpressionParsingUtil>("exprtk");
        std::vector<std::string> varNames;
        std::vector<double> varVals;
        for (const auto& kv : in_vars)
        {
            varNames.emplace_back(kv.first);
            varVals.emplace_back(kv.second);
        }

        if (!parsingUtil->validExpression(in_exprString, varNames))
        {            
            return false;
        }
 
        double evaled = 0.0;
    
        if (!parsingUtil->evaluate(in_exprString, varNames, varVals, evaled))
        {
            return false;
        }

        out_result = evaled;
        return true;
    }
    
}

namespace QuaC { 
std::unique_ptr<HamiltonianTerm> HamiltonianTimeDependentTerm::fromString(const std::string& in_string, const VarsMap& in_vars) 
{
    auto exprStr = removeWhiteSpaces(in_string);
    // Find the special '||' channel separator
    auto separatorPos = exprStr.find("||");
    if (separatorPos == std::string::npos)
    {
        return nullptr;
    }

    const auto channelName = toUpperCase(exprStr.substr(separatorPos + 2));
    // Minimum length: 2
    if (channelName.length() < 2 || (channelName.front() != 'D' && channelName.front() != 'U') || (!isNumberString(channelName.substr(1))))
    {
        return nullptr;
    }
    
    const auto operatorExpression = exprStr.substr(0, separatorPos);

    if (operatorExpression.back() == ')')
    {
        auto splitExprs = UnwrapOpExpresion(operatorExpression);
        if (splitExprs.size() != 2)
        {
            return nullptr;
        }

        auto expr1 = fromString(splitExprs[0] + exprStr.substr(separatorPos), in_vars);
        auto expr2 = fromString(splitExprs[1] + exprStr.substr(separatorPos), in_vars);
        if (expr1 == nullptr || expr2  == nullptr)
        {
            return nullptr;
        }

        std::vector<std::unique_ptr<HamiltonianTerm>> terms;
        terms.emplace_back(std::move(expr1));
        terms.emplace_back(std::move(expr2));

        return std::unique_ptr<HamiltonianTerm>(new HamiltonianSumTerm(std::move(terms)));
    }

    std::vector<QubitOp> operators;
    QubitOp tempOp;
    std::string remainderStr;
    std::string tempStr = operatorExpression;
    while (GetLastOperator(tempStr, tempOp, remainderStr))
    {
        operators.emplace_back(tempOp);
        tempStr = remainderStr;
    }


    double evaled = 0.0;
    if (!TryEvaluateExpression(remainderStr, in_vars, evaled))
    {
        return nullptr;
    }

    // Reverse the vector list since we were parsing operators from the back.
    std::reverse(operators.begin(), operators.end());
    
    return std::make_unique<HamiltonianTimeDependentTerm>(channelName, evaled, operators);
}

std::unique_ptr<HamiltonianTerm> HamiltonianTimeIndependentTerm::fromString(const std::string& in_string, const VarsMap& in_vars)
{
    auto exprStr = removeWhiteSpaces(in_string);

    // Don't process time-dependent terms
    auto separatorPos = exprStr.find("||");
    if (separatorPos != std::string::npos)
    {
        return nullptr;
    }

    if (exprStr.back() == ')')
    {
        auto splitExprs = UnwrapOpExpresion(exprStr);
        if (splitExprs.size() != 2)
        {
            return nullptr;
        }

        auto expr1 = fromString(splitExprs[0], in_vars);
        auto expr2 = fromString(splitExprs[1], in_vars);
        if (expr1 == nullptr || expr2  == nullptr)
        {
            return nullptr;
        }

        std::vector<std::unique_ptr<HamiltonianTerm>> terms;
        terms.emplace_back(std::move(expr1));
        terms.emplace_back(std::move(expr2));

        return std::unique_ptr<HamiltonianTerm>(new HamiltonianSumTerm(std::move(terms)));
    }

    std::vector<QubitOp> operators;

    QubitOp tempOp;
    std::string remainderStr;
    std::string tempStr = exprStr;
    while (GetLastOperator(tempStr, tempOp, remainderStr))
    {
        operators.emplace_back(tempOp);
        tempStr = remainderStr;
    }
 
    double evaled = 0.0;
    if (!TryEvaluateExpression(remainderStr, in_vars, evaled))
    {
        return nullptr;
    }
    
    // Reverse the vector list since we were parsing operators from the back.
    std::reverse(operators.begin(), operators.end());
    
    return std::make_unique<HamiltonianTimeIndependentTerm>(evaled, operators);
}

std::unique_ptr<HamiltonianTerm> HamiltonianSumTerm::fromString(const std::string& in_string, const VarsMap& in_vars)
{
    static const std::string SUM_TERM_PREFIX = "_SUM[";
    auto exprStr = removeWhiteSpaces(in_string);
    if (exprStr.compare(0, SUM_TERM_PREFIX.size(), SUM_TERM_PREFIX) != 0 || exprStr.back() != ']')
    {
        // Not a Sum term
        return nullptr;
    }

    exprStr =  exprStr.substr(SUM_TERM_PREFIX.size());
    exprStr.pop_back();
    
    const auto GetNextCommaSeparatedSubString = [](std::string& io_string) -> std::string {
        auto commaPos = io_string.find(",");
        if (commaPos == std::string::npos)
        {
            return "";
        }
        
        const auto subStr = io_string.substr(0, commaPos);
        io_string = io_string.substr(commaPos + 1);
        return subStr;
    };
    
    
    const auto loopVarName = GetNextCommaSeparatedSubString(exprStr);
    const auto startValStr = GetNextCommaSeparatedSubString(exprStr);
    const auto endValStr = GetNextCommaSeparatedSubString(exprStr);
    const auto loopExpression = exprStr;
    const std::string varFmt = "{" + loopVarName + "}";
    
    if (loopVarName.empty() || startValStr.empty() || endValStr.empty() || loopExpression.empty() || 
            !isNumberString(startValStr) || !isNumberString(endValStr) ||  
            // The expression doesn't contain the loop index var!!!
            loopExpression.find(varFmt) == std::string::npos) 
    {
        return nullptr;
    }

    const int startLoopVal = std::stoi(startValStr);
    const int endLoopVal = std::stoi(endValStr);

    const auto resolveLoopTemplate = [&varFmt](const std::string& in_string, const std::string& in_loopVarName, int in_val) -> std::string {
	    static const std::regex loopValRegex("(\\{.*?\\})");
        std::string result = in_string;
        std::smatch base_match;
        
        std::regex_search(result, base_match, loopValRegex);
        std::unordered_map<std::string, int> loopVarExprs;
        auto searchStart(result.cbegin());
        while(std::regex_search(searchStart, result.cend(), base_match, loopValRegex)) 
        {
            for (size_t i = 0; i < base_match.size(); ++i) 
            {
                loopVarExprs.emplace(base_match[i], 0);
            }
            searchStart = base_match.suffix().first;
        }
        
        // Just one variable: the loop var
        QuaC::VarsMap loopVarMap;
        loopVarMap.emplace(in_loopVarName, in_val);
        
        for (auto& kv : loopVarExprs)
        {
            // Fast path: if it is simply {i}, the assign the value directly, no need to evaluate
            auto loopVarExpr = kv.first;
            if (loopVarExpr == varFmt)
            {
                kv.second = in_val;
            }
            else
            {
                // Something more complex: like {i+1}
                // Remove the curly brackets
                loopVarExpr.erase(0, 1);
                loopVarExpr.pop_back();
                double evalResult = -1;
                if (TryEvaluateExpression(loopVarExpr, loopVarMap, evalResult))
                {
                    kv.second = (int)evalResult; 
                }
            }
        }

        // Get the first occurrence
        const auto replaceTemplateWithValue = [&](const std::string& in_template, int in_value) {
            size_t pos = result.find(in_template);
            const std::string replaceStr = std::to_string(in_value);
            // Repeat till end is reached
            while(pos != std::string::npos)
            {
                // Replace this occurrence 
                result.replace(pos, in_template.size(), replaceStr);
                // Get the next occurrence from the current position
                pos = result.find(in_template, pos + replaceStr.size());
            }
        };
	    
        for (auto& kv : loopVarExprs)
        {
            replaceTemplateWithValue(kv.first, kv.second);
        }

        return result;
    };


    if (startLoopVal > endLoopVal)
    {
        return nullptr;
    }

    const auto resolvedExpression = resolveLoopTemplate(loopExpression, loopVarName, startLoopVal);
    const auto tryTimeIndependent = HamiltonianTimeIndependentTerm::fromString(resolvedExpression, in_vars);
    const auto tryTimeDependent = HamiltonianTimeDependentTerm::fromString(resolvedExpression, in_vars);

    if (tryTimeIndependent == nullptr && tryTimeDependent == nullptr)
    {
        return nullptr;
    }

    const bool isTimeDependent = tryTimeDependent != nullptr;
    const auto parseLoopExpression = [&](const std::string& in_exprStr) ->  std::unique_ptr<HamiltonianTerm> {
        if (isTimeDependent)
        {
            auto result = HamiltonianTimeDependentTerm::fromString(in_exprStr, in_vars);
            return std::unique_ptr<HamiltonianTerm>(result.release()); 
        }
        else
        {
            auto result = HamiltonianTimeIndependentTerm::fromString(in_exprStr, in_vars);
            return std::unique_ptr<HamiltonianTerm>(result.release()); 
        }        
    };

    // Note: IBM uses an inclusive loop index (i.e. the end value is included)
    std::vector<std::unique_ptr<HamiltonianTerm>> loopOps;

    for (int i = startLoopVal; i <= endLoopVal; ++i)
    {
        const auto resolvedExpression = resolveLoopTemplate(loopExpression, loopVarName, i);
        std::unique_ptr<HamiltonianTerm> result = parseLoopExpression(resolvedExpression);
        assert(result != nullptr);
        loopOps.emplace_back(std::move(result));
    }

    return std::make_unique<HamiltonianSumTerm>(std::move(loopOps));
}

void HamiltonianTimeIndependentTerm::apply(IChannelNameResolver* in_channelResolver, FunctorExecutorBase* in_executor)
{
    // This constraint can be lifted if necessary, just add API's to the backend. 
    if (m_operators.size() > 2)
    {
        xacc::error("We only support Hamiltonian terms which are products of maximum two operators.");
    }
    
    if (m_operators.size() == 1)
    {
        const auto op = m_operators.front();          
        std::vector<std::pair<std::string, int>> opList;
        opList.emplace_back(OperatorToString(op.first), op.second);
        
        in_executor->PostFunctorAsync(std::make_unique<AddHamiltonianTerm>(
            m_coefficient,
            opList
        ));      
    }
    else if (m_operators.size() == 2)
    {
        const auto op1 = m_operators[0];
        const auto op2 = m_operators[1];
        std::vector<std::pair<std::string, int>> opList;
        opList.emplace_back(OperatorToString(op1.first), op1.second);
        opList.emplace_back(OperatorToString(op2.first), op2.second);
        in_executor->PostFunctorAsync(std::make_unique<AddHamiltonianTerm>(
            m_coefficient,
            opList
        ));      
    }
}

void HamiltonianTimeDependentTerm::apply(IChannelNameResolver* in_channelResolver, FunctorExecutorBase* in_executor)
{
    // We only support multiplication of up to two operators
    assert(m_operators.size() == 1 || m_operators.size() == 2);
    if (m_operators.size() == 1)
    {
        const auto op = m_operators.front();          
        std::vector<std::pair<std::string, int>> opList;
        opList.emplace_back(OperatorToString(op.first), op.second);
        
        in_executor->PostFunctorAsync(std::make_unique<AddHamiltonianTerm>(
            m_coefficient,
            opList,
            in_channelResolver->GetChannelId(m_channelName)
        ));              
    }
    else if (m_operators.size() == 2)
    {
        const auto op1 = m_operators[0];
        const auto op2 = m_operators[1];
        std::vector<std::pair<std::string, int>> opList;
        opList.emplace_back(OperatorToString(op1.first), op1.second);
        opList.emplace_back(OperatorToString(op2.first), op2.second);
        in_executor->PostFunctorAsync(std::make_unique<AddHamiltonianTerm>(
            m_coefficient,
            opList,
            in_channelResolver->GetChannelId(m_channelName)
        ));
    }    
}

void HamiltonianSumTerm::apply(IChannelNameResolver* in_channelResolver, FunctorExecutorBase* in_executor)
{
    for (auto& term : m_terms)
    {
        term->apply(in_channelResolver, in_executor);
    }
}

void HamiltonianTimeIndependentTerm::collect(std::string& io_staticHstr, std::vector<std::string>& io_ctrlHstr, std::vector<std::string>& io_channelNames)
{
    if (m_operators.size() > 2)
    {
        xacc::error("We only support Hamiltonian terms which are products of maximum two operators.");
    }
    // Don't need to add zero term
    if (std::abs(m_coefficient) < 1e-12)
    {
        return;
    }

    if (m_operators.size() == 1)
    {
        const auto op = m_operators.front();        
        const std::string hStr =  "+ " + std::to_string(m_coefficient.real()) + " " + OperatorToString(op.first) + std::to_string(op.second);
        // Debug:
        std::cout << "add static H string: " << hStr << "\n";
        io_staticHstr.append(hStr);
    }
    else if (m_operators.size() == 2)
    {
        const auto op1 = m_operators[0];
        const auto op2 = m_operators[1];
       
        const std::string hStr =  "+ " + std::to_string(m_coefficient.real()) + " " + 
            OperatorToString(op1.first) + std::to_string(op1.second) +
            OperatorToString(op2.first) + std::to_string(op2.second) ;
        
        // Debug:
        std::cout << "add static H string: " << hStr << "\n";
        io_staticHstr.append(hStr);
    }
}

void HamiltonianTimeDependentTerm::collect(std::string& io_staticHstr, std::vector<std::string>& io_ctrlHstr, std::vector<std::string>& io_channelNames)
{
    // We only support multiplication of up to two operators
    assert(m_operators.size() == 1 || m_operators.size() == 2);
    io_channelNames.emplace_back(m_channelName);
    if (m_operators.size() == 1)
    {
        const auto op = m_operators.front();          
        const std::string hStr =  std::to_string(m_coefficient) + " " + OperatorToString(op.first) + std::to_string(op.second);
        // Debug:
        std::cout << "add control H string: " << hStr << "\n";
        io_ctrlHstr.emplace_back(hStr);      
    }
    else if (m_operators.size() == 2)
    {
        const auto op1 = m_operators[0];
        const auto op2 = m_operators[1];
        
        const std::string hStr =  std::to_string(m_coefficient) + " " + 
            OperatorToString(op1.first) + std::to_string(op1.second) +
            OperatorToString(op2.first) + std::to_string(op2.second) ;
        
        // Debug:
        std::cout << "add control H string: " << hStr << "\n";
        io_ctrlHstr.emplace_back(hStr);      
    }    
}

void HamiltonianSumTerm::collect(std::string& io_staticHstr, std::vector<std::string>& io_ctrlHstr, std::vector<std::string>& io_channelNames)
{
    for (auto& term : m_terms)
    {
        term->collect(io_staticHstr, io_ctrlHstr, io_channelNames);
    }
}

std::unique_ptr<HamiltonianTerm> HamiltonianTimeIndependentTerm::clone()
{
    return std::unique_ptr<HamiltonianTerm>(new HamiltonianTimeIndependentTerm(m_coefficient, m_operators)); 
}

std::unique_ptr<HamiltonianTerm> HamiltonianTimeDependentTerm::clone()
{
    return std::unique_ptr<HamiltonianTerm>(new HamiltonianTimeDependentTerm(m_channelName, m_coefficient, m_operators)); 
}

std::unique_ptr<HamiltonianTerm> HamiltonianSumTerm::clone()
{
    std::vector<std::unique_ptr<HamiltonianTerm>> clones;
    for (auto& term : m_terms)
    {
        clones.emplace_back(term->clone());
    }

    return std::unique_ptr<HamiltonianTerm>(new HamiltonianSumTerm(std::move(clones)));
}

std::unique_ptr<HamiltonianTerm> HamiltonianParsingUtil::tryParse(const std::string& in_expr, const VarsMap& in_vars)
{
    {
        auto trySum = HamiltonianSumTerm::fromString(in_expr, in_vars);
        if (trySum)
        {
            return trySum; 
        }
    }
    
    {
        auto tryTimeDep = HamiltonianTimeDependentTerm::fromString(in_expr, in_vars);
    
        if (tryTimeDep)
        {
            return tryTimeDep;

        }
    }

    {
        auto tryTimeInd = HamiltonianTimeIndependentTerm::fromString(in_expr, in_vars);
    
        if (tryTimeInd)
        {
            return tryTimeInd;
        }
    }
    
    std::cout << "Cannot parse Hamiltonian string " << in_expr << "\n";
    return nullptr;
}

bool HamiltonianParsingUtil::tryParse(const std::string& in_jsonString, std::function<void(HamiltonianTerm&)> in_forEachTermFn)
{
    auto j = nlohmann::json::parse(in_jsonString);
    if (!j.is_object())
    {
        xacc::warning("Hamiltonian JSON must be an object.");
        return false;
    }

    // Hamiltonian strings and vars map
    auto hamStrArray = j["h_str"];
    auto varsArray = j["vars"];
    
    VarsMap vars;
    for (auto varsArrayIter = varsArray.begin(); varsArrayIter != varsArray.end(); ++varsArrayIter)
    {
        vars.emplace(varsArrayIter.key(), varsArrayIter.value().get<double>());
    }

    for (auto hamStrIter = hamStrArray.begin(); hamStrIter != hamStrArray.end(); ++hamStrIter) 
    {
        auto hamStr = (*hamStrIter).get<std::string>();
        std::cout << "Hamiltonian term: " << hamStr << "\n";
        auto parseResult = tryParse(hamStr, vars);
        if (parseResult == nullptr)
        {
            return false;
        }

        in_forEachTermFn(*parseResult);
    }

    return true;
}
}