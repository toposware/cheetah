# Requires to install inflect package first
# Run `pip install inflect`

import re
import inflect
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

##################################
# FIELD AND RING DEFINITIONS
##################################

p_string = "2**64 - 2**32 + 1"
p = Integer(eval(p_string))
Fp = GF(p)
K = Fp['x']
x = K.gen()
poly = x**6 - 7
Fp6 = Fp.extension(poly, 'u')

u = Fp6.gen()

q = 10230194559610033405039867617070259612247645045591847851798073552054039295467
Fq = GF(q)  # scalar field of the curve

p_len_bound = 16 * min(i for i in range(0, 32) if p.nbits() <= i*16)
q_len_bound = 16 * min(i for i in range(0, 32) if q.nbits() <= i*16)

# curve equation constant (y^2 = x^3 + x + B)
B = 8751648879042009540*u ^ 5 + 9152663345541292493*u ^ 4 + 11461655964319546493 * \
    u ^ 3 + 18139339604299112321*u ^ 2 + 11484739320139982777*u + 13239264519837388909

# curve generator coordinates
g_x = 14525701791208974225*u ^ 5 + 9317930114507436616*u ^ 4 + 15366750451009235558 * \
    u ^ 3 + 1487465174323083563*u ^ 2 + 7015637262057020447*u + 7245987309612903154
g_y = 15793083952783722763*u ^ 5 + 5153214154709431636*u ^ 4 + 17741001138250422543 * \
    u ^ 3 + 5112079552571193492 * u ^ 2 + \
    14203600126873428360*u + 5846570486458364631

E = EllipticCurve(Fp6, [1, B])
G = E(g_x, g_y)


def twoadicity(x, base=0, limit=256):
    return max(i for i in range(base, limit) if ((x-1) % (1 << i) == 0))


p_adicity = twoadicity(p)
t_p = (p - 1) // 2 ^ p_adicity
g_p = Fp.multiplicative_generator()
root_p = g_p ^ t_p

t_p6 = (p ^ 6 - 1) // 2 ^ (p_adicity + 1)
g_p6 = Fp6.multiplicative_generator()
root_p6 = g_p6 ^ t_p6

##################################
# FIELD UTILITY FUNCTIONS
##################################


def repr_field_element_bytes(n, output_hex=False, nb_bytes=32):
    assert nb_bytes % 16 == 0
    n = str(hex(Integer(n)))[2:]
    while len(n) < nb_bytes * 2:
        n = "0" + n
    res = []
    for i in range(nb_bytes-1, -1, -1):
        if output_hex:
            res.append("0x" + n[i*2:i*2+2])
        else:
            res.append(Integer("0x" + n[i*2:i*2+2]))
    return res


def repr_arbitrary_element(n, output_hex=False, nb_hex=64):
    assert nb_hex % 16 == 0
    n_hex_len = len(hex(Integer(n)))
    if nb_hex < n_hex_len - 2:
        nb_hex = 4 * min(i for i in range(4, 2 * n_hex_len, 4)
                         if n_hex_len < i*4)
    n = str(hex(Integer(n)))[2:]
    while len(n) < nb_hex:
        n = "0" + n
    num_list = reversed(range(nb_hex//16))

    output = "[\n"
    for i in num_list:
        if output_hex:
            output += f"    0x{n[i*16:i*16+16]},\n"
        else:
            int = Integer("0x" + n[i*16:i*16+16])
            output += f"    {int},\n"
    output += "]"

    return output


def repr_fp(n, output_hex=False, allow_neg=False, no_computations=False):
    assert(n in Fp)
    if n == 0:
        return "Fp::zero()"
    if n == 1:
        return "Fp::one()"
    is_neg = allow_neg and n > (p-1)//2
    output = ""
    if is_neg:
        output += "(&"
        n = p - n
    if no_computations:
        output += "Fp("
    else:
        output += "Fp::new("
    if is_neg:
        output += f"{hex(n) if output_hex else n})).neg()"
    else:
        output += f"{hex(n) if output_hex else n})"

    return output


def repr_fp6(n, output_hex=False, allow_neg=False, no_computations=False):
    assert(n in Fp6)
    if n == 0:
        return "Fp6::zero()"
    if n == 1:
        return "Fp6::one()"
    c = n.polynomial().coefficients(sparse=False)
    while len(c) < 6:
        c.append(0)
    output = "Fp6 {\n"
    output += "c0: {},\n".format(repr_fp(c[0],
                                 output_hex, allow_neg, no_computations))
    output += "c1: {},\n".format(repr_fp(c[1],
                                 output_hex, allow_neg, no_computations))
    output += "c2: {},\n".format(repr_fp(c[2],
                                 output_hex, allow_neg, no_computations))
    output += "c3: {},\n".format(repr_fp(c[3],
                                 output_hex, allow_neg, no_computations))
    output += "c4: {},\n".format(repr_fp(c[4],
                                 output_hex, allow_neg, no_computations))
    output += "c5: {},\n".format(repr_fp(c[5],
                                 output_hex, allow_neg, no_computations))
    output += "}"
    return (output)


def repr_scalar(n, output_hex=True, no_computations=False):
    assert(n in Fq)
    output = ""
    if no_computations:
        output += "Scalar([\n"
        n = Fq(n) * 2 ^ q_len_bound
    else:
        output += "Scalar::new([\n"
    n = str(hex(Integer(n)))[2:]
    while len(n) < q_len_bound//4:
        n = "0" + n
    for i in range((q_len_bound//64) - 1, -1, -1):
        string = "0x" + n[i*16:i*16+16]
        output += f"    {string if output_hex else Integer(string)},\n"

    output += "])"
    return output


def repr_point(p, output_hex=False, allow_neg=False, no_computations=False):
    assert(p in E)
    output = "ProjectivePoint {\n"
    output += "x: {},\n".format(repr_fp6(p.xy()[0],
                                         output_hex, allow_neg, no_computations))
    output += "y: {},\n".format(repr_fp6(p.xy()[1],
                                         output_hex, allow_neg, no_computations))
    output += "z: Fp6::one(),\n"
    output += "}"
    return (output)


def print_field_element_bytes(n, output_hex=False, nb_bytes=32):
    print(repr_field_element_bytes(n, output_hex, nb_bytes))


def print_arbitrary_element(n, output_hex=False, nb_hex=32):
    print(repr_arbitrary_element(n, output_hex, nb_hex))


def print_fp(n, output_hex=False, allow_neg=False, no_computations=False):
    print(repr_fp(n, output_hex, allow_neg, no_computations))


def print_fp6(n, output_hex=False, allow_neg=False, no_computations=False):
    print(repr_fp6(n, output_hex, allow_neg, no_computations))


def print_scalar(n, output_hex=False, no_computations=False):
    print(repr_scalar(n, output_hex, no_computations))


def print_point(p, output_hex=False, allow_neg=False, no_computations=False):
    print(repr_point(p, output_hex, allow_neg, no_computations))


def twoadicity(x):
    return max(i for i in range(0, x.nbits()) if ((x-1) % (1 << i) == 0))


def clean_solution_string(solution):
    """Helper function of solver output"""
    result = []
    inv_list = []
    for member in solution:
        member = str(member)
        member, inv = member.split("/")
        if inv_list == []:
            inv_list.append(inv)
        else:
            assert inv in inv_list

        member = member[7:]
        member = replace_with_self(member)
        member = clean_ints(member)
        result.append(member)

    inv = replace_with_self(inv_list[0][1:-1])
    inv = clean_ints(inv)
    var_names, var_values, inv = remove_duplicates(inv)
    for i, m in enumerate(result):
        _, _, result[i] = remove_duplicates(m)

    result.append(inv)
    result.append(var_names)
    result.append(var_values)
    return result


def replace_with_self(string):
    """Helper function of solver output"""
    string = string.replace("c0", "self.c0")
    string = string.replace("c1", "self.c1")
    string = string.replace("c2", "self.c2")
    string = string.replace("c3", "self.c3")
    string = string.replace("c4", "self.c4")
    string = string.replace("c5", "self.c5")
    string = string.replace("^2", ".square()")

    return string


def clean_ints(string):
    if string[0] == "(" and string[-1] == ")":
        string = string[1:-1]
    l = re.findall(' \d+\*', string) + re.findall('\(\d+\*', string)
    for i in l:
        Int = Fp(int(i[1:-1]))
        is_neg = Int > (p-1)//2
        q = string.split(i)
        string = q[0]
        for j in range(1, len(q)):
            if is_neg:
                if string[-1] == "-":
                    string = string[:-1] + "+"
                    string += f"{i[0]}Fp::new({p-Int})*{q[j]}"
                else:
                    string = string[:-1] + "-"
                    string += f"{i[0]}Fp::new({p-Int})*{q[j]}"
            else:
                string += f"{i[0]}Fp::new({Int})*{q[j]}"
    return string


def remove_duplicates(string):
    l = set(re.findall('Fp::new\(\d+\)', string))
    e = inflect.engine()
    variable_names = []
    variable_values = []
    for i in l:
        integer = re.findall('\d+', i)[0]
        variable = e.number_to_words(integer).replace("-", "_")
        variable_names.append(variable)
        variable_values.append(i)
        string = string.replace(i, variable)

    string = string.replace("self.c0^2", "c0_sq")
    string = string.replace("self.c1^2", "c1_sq")
    string = string.replace("self.c2^2", "c2_sq")
    string = string.replace("self.c3^2", "c3_sq")
    string = string.replace("self.c4^2", "c4_sq")
    string = string.replace("self.c5^2", "c5_sq")

    string = string.replace("self.c0^3", "self.c0 * c0_sq")
    string = string.replace("self.c1^3", "self.c1 * c1_sq")
    string = string.replace("self.c2^3", "self.c2 * c2_sq")
    string = string.replace("self.c3^3", "self.c3 * c3_sq")
    string = string.replace("self.c4^3", "self.c4 * c4_sq")
    string = string.replace("self.c5^3", "self.c5 * c5_sq")

    return variable_names, variable_values, string


##################################
# FIELD CONSTANTS
##################################

def print_fp_constants():
    p_repr = p_string.replace("**", "^")

    output = "\n// ******************************** //\n"
    output += "// ********* FP CONSTANTS ********* //\n"
    output += "// ******************************** //\n"
    output += f"\n// Field modulus = {p_repr}\n"
    output += f"const M: Fp = Fp({p});\n\n"

    output += f"// 2^{p_len_bound} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R: Fp = {repr_fp(Fp(1), False, False, True)};\n\n"

    output += f"// 2^{p_len_bound * 2} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R2: Fp = {repr_fp(Fp(2^p_len_bound), False, False, True)};\n\n"

    output += f"// 2^{p_len_bound * 3} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R3: Fp = {repr_fp(Fp(2^(p_len_bound * 2)), False, False, True)};\n\n"

    output += f"// Multiplicative generator g of order p-1\n"
    output += f"// g = {g_p}\n"
    output += f"//   = {g_p * 2^p_len_bound} in Montgomery form\n"
    output += f"const GENERATOR: Fp = {repr_fp(g_p, False, False, True)};\n\n"

    output += f"// Two-adicity of the field: (p-1) % 2^{p_adicity} = 0\n"
    output += f"pub(crate) const TWO_ADICITY: u32 = {p_adicity};\n\n"

    output += f"// 2^{p_adicity} root of unity = {root_p}\n"
    output += f"//                    = {root_p * 2^p_len_bound} in Montgomery form\n"
    output += f"const TWO_ADIC_ROOT_OF_UNITY: Fp = {repr_fp(root_p, True, False, True)};\n\n"

    output += f"// -M^{{-1}} mod 2^{p_len_bound}; this is used during element multiplication.\n"
    output += f"const U: u64 = {-p^(-1) % 2^p_len_bound};\n\n"

    print(output)


def print_fp6_constants():
    poly = -Fp6.modulus()
    coeffs = poly.coefficients(sparse=False)

    output = "\n// ******************************** //\n"
    output += "// ******** FP6  CONSTANTS ******** //\n"
    output += "// ******************************** //\n\n"

    for i in range(len(coeffs)-1):
        if coeffs[i] != 0:
            output += f"// Necessary for multiplication by v^3\n"
            output += f"const V{i}: {repr_fp(coeffs[i], False, True)};\n\n"

    output += f"// 2^{p_adicity + 1} root of unity = {root_p6}\n"
    output += f"//                    = {root_p6 * 2^p_len_bound} in Montgomery form\n"
    output += f"const TWO_ADIC_ROOT_OF_UNITY_P6: Fp6 = {repr_fp6(root_p6, True, False, True)};\n\n"

    output += f"const TWO_ADICITY_P6: u32 = TWO_ADICITY + 1;\n\n"

    print(output)


def print_scalar_constants():
    s = twoadicity(q)
    t = (q-1) // 2 ^ s
    g = Fq.multiplicative_generator()
    root_of_unity = g ^ t

    output = "\n// ******************************** //\n"
    output += "// ********* FQ CONSTANTS ********* //\n"
    output += "// ******************************** //\n"
    output += f"\n// Field modulus = {q}\n"
    output += f"const M: Scalar = Scalar({repr_arbitrary_element(q, True)});\n\n"

    output += f"// 2^{q_len_bound} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R: Scalar = {repr_scalar(Fq(1), True, True)};\n\n"

    output += f"// 2^{q_len_bound * 2} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R2: Scalar = {repr_scalar(Fq(2^q_len_bound), True, True)};\n\n"

    output += f"// 2^{q_len_bound * 3} mod M; this is used for conversion of elements into Montgomery representation.\n"
    output += f"pub(crate) const R3: Scalar = {repr_scalar(Fq(2^(q_len_bound * 2)), True, True)};\n\n"

    output += f"// Multiplicative generator g of order q-1\n"
    output += f"// g = {g}\n"
    output += f"//   = {hex(g * 2^q_len_bound)} in Montgomery form\n"
    output += f"const GENERATOR: Scalar = {repr_scalar(g, True, True)};\n\n"

    output += f"// Two-adicity of the field: (q-1) % 2^{s} = 0\n"
    output += f"const TWO_ADICITY: u32 = {s};\n\n"

    output += f"// 2^{s} root of unity = {hex(root_of_unity)}\n"
    output += f"//                   = {hex(root_of_unity * 2^q_len_bound)} in Montgomery form\n"
    output += f"const TWO_ADIC_ROOT_OF_UNITY: Scalar = {repr_scalar(root_of_unity, True, True)};\n\n"

    output += f"// -M^{{-1}} mod 2^{p_len_bound}; this is used during element multiplication.\n"
    output += f"const U: u64 = {Fq(-q^(-1) % 2^p_len_bound)};\n\n"

    print(output)

##################################
# FIELD INVERSION FUNCTIONS
##################################


def find_inverse_formula_fp6():
    """Non-optimized formula for inversion in Fp6"""
    # solving in \mathbf{Z} is enough
    c0, c1, c2, c3, c4, c5, c0p, c1p, c2p, c3p, c4p, c5p = var(
        "c0, c1, c2, c3, c4, c5, c0p, c1p, c2p, c3p, c4p, c5p", domain=ZZ)
    t0, t1, t2, t3, t4 = var("t0, t1, t2, t3, t4", domain=ZZ)
    poly = Fp6.modulus()
    (u0, u1, u2, u3, u4, u5, _u6) = (-poly).coefficients(sparse=False)
    (u0, u1, u2, u3, u4, u5) = int(u0), int(
        u1), int(u2), int(u3), int(u4), int(u5)

    solution = solve([c1 * c5p + c2 * c4p + c3 * c3p + c4 * c2p + c5 * c1p == t0,
                      c2 * c5p + c3 * c4p + c4 * c3p + c5 * c2p == t1,
                      c3 * c5p + c4 * c4p + c5 * c3p == t2,
                      c4 * c5p + c5 * c4p == t3,
                      c5 * c5p == t4,
                      c0 * c0p + t0 == 1,
                      c0 * c1p + c1 * c0p + t0 + t1 == 0,
                      c0 * c2p + c1 * c1p + c2 * c0p + t1 + t2 == 0,
                      c0 * c3p + c1 * c2p + c2 * c1p + c3 * c0p + t2 + t3 == 0,
                      c0 * c4p + c1 * c3p + c2 * c2p + c3 * c1p + c4 * c0p + t3 + t4 == 0,
                      c0 * c5p + c1 * c4p + c2 * c3p + c3 * c2p + c4 * c1p + c5 * c0p + t4 == 0, ], c0p, c1p, c2p, c3p, c4p, c5p)[0]
    c0p, c1p, c2p, c3p, c4p, c5p, inv, var_names, var_values = clean_solution_string(
        solution)

    output = "\n// ******************************** //\n"
    output += "// ******** FP6  INVERSION ******** //\n"
    output += "// ******************************** //\n\n"

    for i, v in enumerate(var_names):
        output += f"let {v} = {var_values[i]};\n"

    for c in ["c0", "c1", "c2", "c3", "c4", "c5"]:
        if c + "_sq" in inv:
            output += f"\nlet {c}_sq = self.{c}.square();"

    output += f"\n\nlet inv = {inv};\n\n"
    output += f"let c0 = {c0p};\n"
    output += f"let c1 = {c1p};\n"
    output += f"let c2 = {c2p};\n\n"
    output += f"let c3 = {c3p};\n\n"
    output += f"let c4 = {c4p};\n\n"
    output += f"let c5 = {c5p};\n\n"
    output += f"inv.invert().map( | t | Fp6 {{\n    c0: c0 * t,\n    c1: c1 * t,\n    c2: c2 * t,\n    c3: c3 * t,\n    c4: c4 * t,\n    c5: c5 * t,\n}})\n"

    output = output.replace("Fp::", "Fp2::")
    print(output)


##################################
# FIELD MULTIPLICATION FUNCTIONS
##################################

def print_multiplication_formula_fp6():
    """Generic non-optimized formula for multiplication in Fp6"""
    poly = Fp6.modulus()
    (u0, u1, u2, u3, u4, u5, _u6) = (-poly).coefficients(sparse=False)
    (u0, u1, u2, u3, u4, u5) = int(u0), int(
        u1), int(u2), int(u3), int(u4), int(u5)

    output = "\n// ******************************** //\n"
    output += "// ****** FP6 MULTIPLICATION ****** //\n"
    output += "// ******************************** //\n\n"

    output += "let aa = (&self.c0).mul(&other.c0);\n"
    output += "let ab = (&self.c0).mul(&other.c1);\n"
    output += "let ac = (&self.c0).mul(&other.c2);\n"
    output += "let ad = (&self.c0).mul(&other.c3);\n"
    output += "let ae = (&self.c0).mul(&other.c4);\n"
    output += "let af = (&self.c0).mul(&other.c5);\n\n"

    output += "let ba = (&self.c1).mul(&other.c0);\n"
    output += "let bb = (&self.c1).mul(&other.c1);\n"
    output += "let bc = (&self.c1).mul(&other.c2);\n"
    output += "let bd = (&self.c1).mul(&other.c3);\n"
    output += "let be = (&self.c1).mul(&other.c4);\n"
    output += "let bf = (&self.c1).mul(&other.c5);\n\n"

    output += "let ca = (&self.c2).mul(&other.c0);\n"
    output += "let cb = (&self.c2).mul(&other.c1);\n"
    output += "let cc = (&self.c2).mul(&other.c2);\n"
    output += "let cd = (&self.c2).mul(&other.c3);\n"
    output += "let ce = (&self.c2).mul(&other.c4);\n"
    output += "let cf = (&self.c2).mul(&other.c5);\n\n"

    output += "let da = (&self.c3).mul(&other.c0);\n"
    output += "let db = (&self.c3).mul(&other.c1);\n"
    output += "let dc = (&self.c3).mul(&other.c2);\n"
    output += "let dd = (&self.c3).mul(&other.c3);\n"
    output += "let de = (&self.c3).mul(&other.c4);\n"
    output += "let df = (&self.c3).mul(&other.c5);\n\n"

    output += "let ea = (&self.c4).mul(&other.c0);\n"
    output += "let eb = (&self.c4).mul(&other.c1);\n"
    output += "let ec = (&self.c4).mul(&other.c2);\n"
    output += "let ed = (&self.c4).mul(&other.c3);\n"
    output += "let ee = (&self.c4).mul(&other.c4);\n"
    output += "let ef = (&self.c4).mul(&other.c5);\n\n"

    output += "let fa = (&self.c5).mul(&other.c0);\n"
    output += "let fb = (&self.c5).mul(&other.c1);\n"
    output += "let fc = (&self.c5).mul(&other.c2);\n"
    output += "let fd = (&self.c5).mul(&other.c3);\n"
    output += "let fe = (&self.c5).mul(&other.c4);\n"
    output += "let ff = (&self.c5).mul(&other.c5);\n\n"

    # u0 cannot be zero
    output += "let t0 = (&bf).add(&fb);\n"
    output += "let t0 = (&t0).add(&ce);\n"
    output += "let t0 = (&t0).add(&ec);\n"
    output += "let t0 = (&t0).add(&dd);\n"
    output += "let c0 = (&t0).add(&aa);\n\n"

    output += "let c1 = (&ab).add(&ba);"
    output += "let c1 = (&c1).add(&t0);"
    output += "let t1 = (&cf).add(&fc);"
    output += "let t1 = (&t1).add(&de);"
    output += "let t1 = (&t1).add(&ed);"
    output += "let c1 = (&c1).add(&t1);"

    output += "let c2 = (&ac).add(&ca);"
    output += "let c2 = (&c2).add(&bb);"
    output += "let c2 = (&c2).add(&t1);"
    output += "let t2 = (&df).add(&fd);"
    output += "let t2 = (&t2).add(&ee);"
    output += "let c2 = (&c2).add(&t2);"

    output += "let c3 = (&ad).add(&da);"
    output += "let c3 = (&c3).add(&bc);"
    output += "let c3 = (&c3).add(&cb);"
    output += "let c3 = (&c3).add(&t2);"
    output += "let t3 = (&ef).add(&fe);"
    output += "let c3 = (&c3).add(&t3);"

    output += "let c4 = (&ae).add(&ea);"
    output += "let c4 = (&c4).add(&bd);"
    output += "let c4 = (&c4).add(&db);"
    output += "let c4 = (&c4).add(&cc);"
    output += "let c4 = (&c4).add(&t3);"
    output += "let c4 = (&c4).add(&ff);"

    output += "let c5 = (&af).add(&fa);"
    output += "let c5 = (&c5).add(&be);"
    output += "let c5 = (&c5).add(&eb);"
    output += "let c5 = (&c5).add(&cd);"
    output += "let c5 = (&c5).add(&dc);"
    output += "let c5 = (&c5).add(&ff);"

    output += f"Fp6 {{ c0, c1, c2, c3, c4, c5 }}\n"
    print(output)


##################################
# FIELD SQUARING FUNCTIONS
##################################

def print_squaring_formula_fp6():
    """Generic non-optimized formula for squaring in Fp6"""
    poly = Fp6.modulus()
    (u0, u1, u2, u3, u4, u5, _u6) = (-poly).coefficients(sparse=False)
    (u0, u1, u2, u3, u4, u5) = int(u0), int(
        u1), int(u2), int(u3), int(u4), int(u5)

    output = "\n// ******************************** //\n"
    output += "// ********* FP6 SQUARING ********* //\n"
    output += "// ******************************** //\n\n"

    output += "let aa = (&self.c0).mul(&self.c0);\n"
    output += "let ab = (&self.c0).mul(&self.c1);\n"
    output += "let ac = (&self.c0).mul(&self.c2);\n\n"
    output += "let bb = (&self.c1).mul(&self.c1);\n"
    output += "let bc = (&self.c1).mul(&self.c2);\n\n"
    output += "let cc = (&self.c2).mul(&self.c2);\n\n"

    # u0 cannot be zero
    output += "let c0 = (&bc).double();\n"
    if u2 != 0:
        output += "let t0 = (&V2).mul(&cc);\n"
        output += "let c0 = (&c0).add(&t0);\n"
    output += "let c0 = (&c0).mul(&V0);\n"
    output += "let c0 = (&c0).add(&aa);\n\n"

    if u1 != 0:
        output += "let c1 = (&bc).double();\n"
        output += "let c1 = (&c1).mul(&V1);\n"
        if u2 != 0:
            output += "let t1 = (&V1).mul(&V2);\n"
            output += "let t1 = (&t1).add(&V0);\n"
            output += "let t1 = (&t1).mul(&cc);\n"
        else:
            output += "let t1 = (&cc).mul(&V0);\n"
        output += "let c1 = (&c1).add(&t1);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ab);\n\n"
    else:
        output += "let c1 = (&V0).mul(&cc);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ab);\n\n"

    if u2 != 0:
        output += "let c2 = (&bc).double();\n"
        output += "let c2 = (&c2).mul(&V2);\n"
        output += "let t2 = (&V2).square();\n"
        if u1 != 0:
            output += "let t2 = (&t2).add(&V1);\n"
        output += "let t2 = (&t2).mul(&cc);\n"
        output += "let c2 = (&c2).add(&t2);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&bb);\n\n"
    else:
        output += "let c2 = (&V1).mul(&cc);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&bb);\n\n"

    output += f"Fp6 {{ c0, c1, c2 }}\n"
    print(output)


def print_basepoint_table():
    point = G
    L = []
    for i in range(32):
        M = []
        for j in range(8):
            M.append(point)
            point += G
        L.append(M)
        point = 256 * point

    output = "BasePointTable([\n"
    for table in L:
        output += "\tLookupTable([\n"
        for t in table:
            output += "\t\t"
            output += repr_point(t, True, False, True)
            output += ",\n"
        output += "]),\n"
    output += "]);"
    print(output)


# print_fp_constants()

# print_scalar_constants()

# print_fp6_constants()
# print_multiplication_formula_fp6()
# print_squaring_formula_fp6()
# find_inverse_formula_fp6()

print_basepoint_table()
