# Requires to install inflect package first
# Run `pip install inflect`

import re
import inflect
import warnings
warnings.simplefilter("ignore", category=DeprecationWarning)

##################################
# FIELD AND RING DEFINITIONS
##################################

p_string = "2**62 + 2**56 + 2**55 + 1"
p = Integer(eval(p_string))
Fp = GF(p)
K = Fp['x']
x = K.gen()
poly2 = x**2 - 2*x - 2
Fp2 = Fp.extension(poly2, 'u')
K2 = Fp2['x']
x = K2.gen()
poly3 = x**3 + x + 1
Fp6 = Fp2.extension(poly3, 'v')

u = Fp2.gen()
v = Fp6.gen()

q = 17278877126736494933592566161653303514319447234579276854188469089485337225893
Fq = GF(q)  # scalar field of the curve

p_len_bound = 16 * min(i for i in range(0, 32) if p.nbits() < i*16)
q_len_bound = 16 * min(i for i in range(0, 32) if q.nbits() < i*16)

# curve equation constant (y^2 = x^3 + x + B)
B = (1200866201009650596*u + 1935817186716799185)*v ^ 2 + (3999205700308519553 *
                                                           u + 3518137720867787056)*v + 2508413708960025374*u + 1526905369741321712

# curve generator coordinates
g_x = (288076929228681448*u + 2633256936270674947)*v ^ 2 + (1056103921720638754 *
                                                            u + 3052857668015466949)*v + 4508025770867562887*u + 2398517019392108645
g_y = (3289504647774244396*u + 4227116334258416103)*v ^ 2 + (3024200307602630234 *
                                                             u + 1961556908722893436)*v + 1225290585625954719*u + 3894155704139868264

E = EllipticCurve(Fp6, [1, B])


def make_finite_field(k):
    r""" Return the finite field isomorphic to this field.

    INPUT:

    - ``k`` -- a finite field

    OUTPUT: a tuple `(k_1,\phi,\xi)` where `k_1` is a 'true' finite field,
    `\phi` is an isomorphism from `k` to `k_1` and `\xi` is an isomorphism
    from `k_1` to `k`.

    This function is useful when `k` is constructed as a tower of extensions
    with a finite field as a base field.

    Adapted from https://github.com/MCLF/mclf/issues/103.

    """

    assert k.is_field()
    assert k.is_finite()
    # TODO: partially solved sage9.4 issue but still failing for higher extensions (wrong isomorphic field)
    if k.base_ring().is_prime_field():
        return k, k.hom(k.gen(), k), k.hom(k.gen(), k)
    else:
        k0 = k.base_field()
        G = k.modulus()
        assert G.parent().base_ring() is k0
        k0_new, phi0, _ = make_finite_field(k0)
        G_new = G.map_coefficients(phi0, k0_new)
        k_new = k0_new.extension(G_new.degree())

        alpha = G_new.roots(k_new)[0][0]
        Pk0 = k.cover_ring()
        Pk0_new = k0_new[Pk0.variable_name()]
        psi1 = Pk0.hom(phi0, Pk0_new)
        psi2 = Pk0_new.hom(alpha, k_new)
        psi = psi1.post_compose(psi2)
        # psi: Pk0 --> k_new
        phi = k.hom(Pk0.gen(), Pk0, check=False)
        phi = phi.post_compose(psi)

        k_inv = k0.base_ring()
        phi0_inv = k_inv.hom(k_inv.gen(), k_inv)
        G_new_inv = k_new.modulus().map_coefficients(phi0_inv, k0_new)
        alpha_inv = G_new_inv.roots(k)[0][0]
        phi_inv = k_new.hom(alpha_inv, k)

        return k_new, phi, phi_inv


def twoadicity(x, base=0, limit=256):
    return max(i for i in range(base, limit) if ((x-1) % (1 << i) == 0))


# Direct extension rather than towered one
# and associated morphisms to and from Fp6
Fp6_direct, phi, psi = make_finite_field(Fp6)

p_adicity = twoadicity(p)
t_p = (p - 1) // 2 ^ p_adicity
g_p = Fp.multiplicative_generator()
root_p = g_p ^ t_p

t_p2 = (p ^ 2 - 1) // 2 ^ (p_adicity + 1)
g_p2 = Fp2.multiplicative_generator()
root_p2 = g_p2 ^ t_p2

t_p6 = (p ^ 6 - 1) // 2 ^ (p_adicity + 1)
g_p6 = Fp6_direct.multiplicative_generator()
root_p6_direct = g_p6 ^ t_p6
root_p6 = psi(root_p6_direct)

##################################
# FIELD UTILITY FUNCTIONS
##################################


def repr_arbitrary_element(n, output_hex=True, nb_hex=64):
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
            output += f"    {Integer(n[i*16:i*16+16])},\n"
    output += "]"

    return output


def repr_fp(n, output_hex=False, allow_neg=False, no_computations=False):
    assert(n in Fp)
    is_neg = allow_neg and n > (p-1)//2
    output = ""
    if is_neg:
        output += "(&"
        n = p - n
    if no_computations:
        output += "Fp("
        n = Fp(n) * 2 ^ p_len_bound
    else:
        output += "Fp::new("
    if is_neg:
        output += f"{hex(n) if output_hex else n})).neg()"
    else:
        output += f"{hex(n) if output_hex else n})"

    return output


def repr_fp2(n, output_hex=False, allow_neg=False, no_computations=False):
    assert(n in Fp2)
    c = (n + Fp2.gen()).polynomial().coefficients(sparse=False)
    c[-1] -= 1
    output = "Fp2 {\n"
    output += "c0: {},\n".format(repr_fp(c[0],
                                 output_hex, allow_neg, no_computations))
    output += "c1: {},\n".format(repr_fp(c[1],
                                 output_hex, allow_neg, no_computations))
    output += "}"
    return (output)


def repr_fp6(n, output_hex=False, allow_neg=False, no_computations=False):
    assert(n in Fp6)
    c = n.list()
    output = "Fp6 {\n"
    output += "c0: {},\n".format(repr_fp2(c[0],
                                 output_hex, allow_neg, no_computations))
    output += "c1: {},\n".format(repr_fp2(c[1],
                                 output_hex, allow_neg, no_computations))
    output += "c2: {},\n".format(repr_fp2(c[2],
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


def print_fp(n, output_hex=False, allow_neg=False, no_computations=False):
    print(repr_fp(n, output_hex, allow_neg, no_computations))


def print_fp2(n, output_hex=False, allow_neg=False, no_computations=False):
    print(repr_fp2(n, output_hex, allow_neg, no_computations))


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

    string = string.replace("self.c0^3", "self.c0 * c0_sq")
    string = string.replace("self.c1^3", "self.c1 * c1_sq")
    string = string.replace("self.c2^3", "self.c2 * c2_sq")

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


def print_fp2_constants():
    poly = -Fp2.modulus()
    coeffs = poly.coefficients(sparse=False)

    output = "\n// ******************************** //\n"
    output += "// ******** FP2  CONSTANTS ******** //\n"
    output += "// ******************************** //\n\n"

    for i in range(len(coeffs)-1):
        if coeffs[i] != 0:
            output += f"// Necessary for multiplication by u^2\n"
            output += f"const U{i}: {repr_fp(coeffs[i], False, True)};\n\n"

    output += f"// 2^{p_adicity + 1} root of unity = {root_p2}\n"
    output += f"//                    = {root_p2 * 2^p_len_bound} in Montgomery form\n"
    output += f"const TWO_ADIC_ROOT_OF_UNITY_P2: Fp2 = {repr_fp2(root_p2, True, False, True)};\n\n"

    output += f"const TWO_ADICITY_P2: u32 = TWO_ADICITY + 1;\n\n"

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
            output += f"const V{i}: {repr_fp2(coeffs[i], False, True)};\n\n"

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


def find_inverse_formula_fp2():
    """Non-optimized formula for inversion in Fp2"""
    # solving in \mathbf{Z} is enough
    c0, c1, c0p, c1p = var("c0, c1, c0p, c1p", domain=ZZ)
    poly = Fp2.modulus()
    (u0, u1, _u2) = (-poly).coefficients(sparse=False)
    (u0, u1) = int(u0), int(u1)

    solution = solve([c1 * c1p * u0 + c0 * c0p == 1,
                      c0 * c1p + c1 * c0p + c1 * c1p * u1 == 0
                      ], c0p, c1p)[0]
    c0p, c1p, inv, var_names, var_values = clean_solution_string(solution)

    output = "\n// ******************************** //\n"
    output += "// ******** FP2  INVERSION ******** //\n"
    output += "// ******************************** //\n\n"

    for i, v in enumerate(var_names):
        output += f"let {v} = {var_values[i]};\n"

    for c in ["c0", "c1"]:
        if c + "_sq" in inv:
            output += f"\nlet {c}_sq = self.{c}.square();"

    output += f"\n\nlet inv = {inv};\n\n"
    output += f"inv.invert().map( | t | Fp2 {{\n    c0: ({c0p}) * t,\n    c1: ({c1p}) * t,\n}})\n"
    print(output)


def find_inverse_formula_fp6():
    """Non-optimized formula for inversion in Fp6"""
    # solving in \mathbf{Z} is enough
    c0, c1, c2, c0p, c1p, c2p = var("c0, c1, c2, c0p, c1p, c2p", domain=ZZ)
    poly = Fp6.modulus()
    (u0, u1, u2, _u3) = (-poly).coefficients(sparse=False)
    (u0, u1, u2) = int(u0), int(u1), int(u2)

    solution = solve([c0 * c0p + u0 * c1*c2p + u0 * c2*c1p + u0*u2*c2*c2p == 1,
                      c0*c1p + c1*c0p + u1*c1*c2p + u1*c2 *
                      c1p + (u0 + u1*u2) * c2*c2p == 0,
                      c0*c2p + c1*c1p + u2*c1*c2p + c2*c0p +
                      u2*c2*c1p + (u1 + u2 ^ 2) * c2*c2p == 0
                      ], c0p, c1p, c2p)[0]
    c0p, c1p, c2p, inv, var_names, var_values = clean_solution_string(solution)

    output = "\n// ******************************** //\n"
    output += "// ******** FP6  INVERSION ******** //\n"
    output += "// ******************************** //\n\n"

    for i, v in enumerate(var_names):
        output += f"let {v} = {var_values[i]};\n"

    for c in ["c0", "c1", "c2"]:
        if c + "_sq" in inv:
            output += f"\nlet {c}_sq = self.{c}.square();"

    output += f"\n\nlet inv = {inv};\n\n"
    output += f"let c0 = {c0p};\n"
    output += f"let c1 = {c1p};\n"
    output += f"let c2 = {c2p};\n\n"
    output += f"inv.invert().map( | t | Fp6 {{\n    c0: c0 * t,\n    c1: c1 * t,\n    c2: c2 * t,\n}})\n"

    output = output.replace("Fp::", "Fp2::")
    print(output)


##################################
# FIELD MULTIPLICATION FUNCTIONS
##################################

def print_multiplication_formula_fp2():
    """Generic non-optimized formula for multiplication in Fp2"""
    poly = Fp2.modulus()
    (_u0, u1, _u2) = (-poly).coefficients(sparse=False)
    u1 = int(u1)

    output = "\n// ******************************** //\n"
    output += "// ****** FP2 MULTIPLICATION ****** //\n"
    output += "// ******************************** //\n\n"

    output += "let aa = (&self.c0).mul(&other.c0);\n"
    output += "let ab = (&self.c0).mul(&other.c1);\n\n"
    output += "let ba = (&self.c1).mul(&other.c0);\n"
    output += "let bb = (&self.c1).mul(&other.c1);\n\n"

    # u0 cannot be zero
    output += "let c0 = (&U0).mul(&bb);\n"
    output += "let c0 = (&c0).add(&aa);\n\n"

    if u1 != 0:
        output += "let c1 = (&U1).mul(&bb);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ba);\n"
    else:
        output += "let c1 = (&ab).add(&ba);\n"

    output += f"\nFp2 {{ c0, c1 }}\n"
    print(output)


def print_multiplication_formula_fp6():
    """Generic non-optimized formula for multiplication in Fp6"""
    poly = Fp6.modulus()
    (_u0, u1, u2, _u3) = (-poly).coefficients(sparse=False)
    (u1, u2) = int(u1), int(u2)

    output = "\n// ******************************** //\n"
    output += "// ****** FP6 MULTIPLICATION ****** //\n"
    output += "// ******************************** //\n\n"

    output += "let aa = (&self.c0).mul(&other.c0);\n"
    output += "let ab = (&self.c0).mul(&other.c1);\n"
    output += "let ac = (&self.c0).mul(&other.c2);\n\n"
    output += "let ba = (&self.c1).mul(&other.c0);\n"
    output += "let bb = (&self.c1).mul(&other.c1);\n"
    output += "let bc = (&self.c1).mul(&other.c2);\n\n"
    output += "let ca = (&self.c2).mul(&other.c0);\n"
    output += "let cb = (&self.c2).mul(&other.c1);\n"
    output += "let cc = (&self.c2).mul(&other.c2);\n\n"

    # u0 cannot be zero
    output += "let c0 = (&bc).add(&cb);\n"
    if u2 != 0:
        output += "let t0 = (&V2).mul(&cc);\n"
        output += "let c0 = (&c0).add(&t0);\n"
    output += "let c0 = (&c0).mul(&V0);\n"
    output += "let c0 = (&c0).add(&aa);\n\n"

    if u1 != 0:
        output += "let c1 = (&bc).add(&cb);\n"
        output += "let c1 = (&c1).mul(&V1);\n"
        if u2 != 0:
            output += "let t1 = (&V1).mul(&V2);\n"
            output += "let t1 = (&t1).add(&V0);\n"
            output += "let t1 = (&t1).mul(&cc);\n"
        else:
            output += "let t1 = (&cc).mul(&V0);\n"
        output += "let c1 = (&c1).add(&t1);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ba);\n\n"
    else:
        output += "let c1 = (&V0).mul(&cc);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ba);\n\n"

    if u2 != 0:
        output += "let c2 = (&bc).add(&cb);\n"
        output += "let c2 = (&c2).mul(&V2);\n"
        output += "let t2 = (&V2).square();\n"
        if u1 != 0:
            output += "let t2 = (&t2).add(&V1);\n"
        output += "let t2 = (&t2).mul(&cc);\n"
        output += "let c2 = (&c2).add(&t2);\n"
        output += "let c2 = (&c2).add(&ca);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&bb);\n\n"
    else:
        output += "let c2 = (&V1).mul(&cc);\n"
        output += "let c2 = (&c2).add(&ca);\n"
        output += "let c2 = (&c2).add(&ac);\n"
        output += "let c2 = (&c2).add(&bb);\n\n"

    output += f"Fp6 {{ c0, c1, c2 }}\n"
    print(output)


##################################
# FIELD SQUARING FUNCTIONS
##################################

def print_squaring_formula_fp2():
    """Generic non-optimized formula for squaring in Fp2"""
    poly = Fp2.modulus()
    (_u0, u1, _u2) = (-poly).coefficients(sparse=False)
    u1 = int(u1)

    output = "\n// ******************************** //\n"
    output += "// ********* FP2 SQUARING ********* //\n"
    output += "// ******************************** //\n\n"

    output += "let aa = &self.c0.square();\n"
    output += "let ab = (&self.c0).mul(&self.c1);\n\n"
    output += "let bb = &self.c1.square();\n\n"

    # u0 cannot be zero
    output += "let c0 = (&U0).mul(&bb);\n"
    output += "let c0 = (&c0).add(&aa);\n\n"

    if u1 != 0:
        output += "let c1 = (&U1).mul(&bb);\n"
        output += "let c1 = (&c1).add(&ab);\n"
        output += "let c1 = (&c1).add(&ab);\n"
    else:
        output += "let c1 = (&ab).double();\n"

    output += f"\nFp2 {{ c0, c1 }}\n"
    print(output)


def print_squaring_formula_fp6():
    """Generic non-optimized formula for squaring in Fp6"""
    poly = Fp6.modulus()
    (_u0, u1, u2, _u3) = (-poly).coefficients(sparse=False)
    (u1, u2) = int(u1), int(u2)

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


print_fp_constants()

print_scalar_constants()

print_fp2_constants()
print_multiplication_formula_fp2()
print_squaring_formula_fp2()
find_inverse_formula_fp2()

print_fp6_constants()
print_multiplication_formula_fp6()
print_squaring_formula_fp6()
find_inverse_formula_fp6()
